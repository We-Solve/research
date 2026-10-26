#include <algorithm>
#include <array>
#include <cassert>
#include <chrono>
#include <cstdint>
#include <fstream>
#include <functional>
#include <iostream>
#include <map>
#include <numeric>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

// Independent exact audit of the reduced k=5 domain in We-Solve/research #17.
// Build: g++ -std=c++17 -O2 -Wall -Wextra -pedantic audit.cpp -o audit
// Run from the directory where the output JSON/JSONL files should be saved.
// No external libraries, floating-point feasibility decisions, or graph atlas.

using Mask = std::uint32_t;
constexpr int K = 5;

int pc(Mask x) { return __builtin_popcount(x); }

struct Graph {
    int n;
    Mask mask;
    std::vector<std::pair<int,int>> edges;
    std::array<Mask,7> adj{};
    std::array<int,7> degree{};
    explicit Graph(int vertices, Mask bits) : n(vertices), mask(bits) {
        int bit = 0;
        for (int u=0; u<n; ++u) for (int v=u+1; v<n; ++v, ++bit) {
            if (mask & (Mask(1)<<bit)) {
                edges.emplace_back(u,v);
                adj[u] |= Mask(1)<<v;
                adj[v] |= Mask(1)<<u;
                ++degree[u]; ++degree[v];
            }
        }
        assert(bit<=21);
        assert((mask >> bit)==0);
    }
};

// Relabel vertices into nondecreasing degree order. Among all permutations
// within equal-degree blocks, take the exact minimum adjacency mask. Every
// graph isomorphism preserves degree; hence this is a complete canonical form.
Mask canonical(const Graph& g) {
    if (g.n<=1) return 0;
    std::vector<std::vector<int>> blocks;
    for (int d=0; d<g.n; ++d) {
        std::vector<int> b;
        for (int v=0; v<g.n; ++v) if (g.degree[v]==d) b.push_back(v);
        if (!b.empty()) blocks.push_back(b);
    }
    std::array<int,7> order{};
    Mask best=~Mask(0);
    std::function<void(std::size_t,int)> visit = [&](std::size_t block, int offset) {
        if (block==blocks.size()) {
            Mask m=0; int bit=0;
            for (int u=0; u<g.n; ++u) for (int v=u+1; v<g.n; ++v, ++bit)
                if (g.adj[order[u]] & (Mask(1)<<order[v])) m |= Mask(1)<<bit;
            best=std::min(best,m);
            return;
        }
        auto b=blocks[block];
        do {
            std::copy(b.begin(),b.end(),order.begin()+offset);
            visit(block+1,offset+static_cast<int>(b.size()));
        } while (std::next_permutation(b.begin(),b.end()));
    };
    visit(0,0);
    return best;
}

Mask extend(const Graph& old, Mask neighbors) {
    const int n=old.n+1;
    Mask m=0; int bit=0;
    for (int u=0; u<n; ++u) for (int v=u+1; v<n; ++v, ++bit) {
        const bool present = v==old.n ? (neighbors & (Mask(1)<<u))
                                             : (old.adj[u] & (Mask(1)<<v));
        if (present) m |= Mask(1)<<bit;
    }
    return m;
}

bool connected_after_removing(const Graph& g, int removed) {
    Mask allowed=(Mask(1)<<g.n)-1;
    if (removed>=0) allowed &= ~(Mask(1)<<removed);
    if (!allowed) return true;
    Mask reached=allowed & (~allowed+1);
    for (;;) {
        Mask next=reached;
        for (int v=0; v<g.n; ++v) if (reached & (Mask(1)<<v))
            next |= g.adj[v] & allowed;
        if (next==reached) break;
        reached=next;
    }
    return reached==allowed;
}

bool core_filter(const Graph& g) {
    const int m=static_cast<int>(g.edges.size());
    if (g.n<5 || g.n>7) return false;
    for (int v=0; v<g.n; ++v) if (g.degree[v]<3) return false;
    if (2*m>K*(g.n-1) || K*(g.n-1)>3*m) return false;
    if (!connected_after_removing(g,-1)) return false;
    for (int v=0; v<g.n; ++v) if (!connected_after_removing(g,v)) return false;
    return true;
}

struct Subset {
    Mask vertices;
    Mask edges;
    int rank_bound;
};

std::vector<Subset> subsets(const Graph& g) {
    std::vector<Subset> result;
    for (Mask s=1; s<(Mask(1)<<g.n)-1; ++s) {
        int size=pc(s);
        if (size<2) continue;
        Mask edges=0;
        for (std::size_t i=0; i<g.edges.size(); ++i) {
            auto [u,v]=g.edges[i];
            if ((s & (Mask(1)<<u)) && (s & (Mask(1)<<v))) edges |= Mask(1)<<i;
        }
        result.push_back({s,edges,size-1});
    }
    return result;
}

template<class F> void combinations(int m, int q, F callback) {
    if (q<0 || q>m) return;
    std::function<void(int,int,Mask)> visit = [&](int start, int left, Mask mask) {
        if (left==0) { callback(mask); return; }
        for (int i=start; i<=m-left; ++i) visit(i+1,left-1,mask|(Mask(1)<<i));
    };
    visit(0,q,0);
}

bool is_tree(const Graph& g, Mask tree) {
    if (pc(tree)!=g.n-1 || (tree >> g.edges.size())!=0) return false;
    std::array<int,7> parent;
    std::iota(parent.begin(),parent.end(),0);
    auto root=[&](int x) { while (parent[x]!=x) x=parent[x]; return x; };
    for (std::size_t i=0; i<g.edges.size(); ++i) if (tree & (Mask(1)<<i)) {
        auto [u,v]=g.edges[i];
        u=root(u); v=root(v);
        if (u==v) return false;
        parent[u]=v;
    }
    for (int v=1; v<g.n; ++v) if (root(v)!=root(0)) return false;
    return true;
}

struct Tree {
    Mask mask;
    std::vector<int> denominators;
};

std::vector<Tree> trees(const Graph& g, const std::vector<Subset>& sets) {
    std::vector<Tree> result;
    combinations(static_cast<int>(g.edges.size()),g.n-1,[&](Mask mask) {
        if (!is_tree(g,mask)) return;
        Tree t{mask,{}};
        for (const auto& s:sets) {
            int d=s.rank_bound-pc(mask&s.edges);
            assert(d>=0);
            t.denominators.push_back(d);
        }
        result.push_back(std::move(t));
    });
    assert(!result.empty());
    return result;
}

struct Rational {
    int numerator;
    int denominator;
    void minimize(int a, int b) {
        assert(a>=0 && b>0);
        if (a*denominator < numerator*b) { numerator=a; denominator=b; }
    }
};

// This direct validation intentionally does not use the optimum-ratio formula.
// It checks tree validity, every residual coordinate and subset inequality,
// then a concrete blocking constraint that forbids every larger coefficient.
// Return (kind,index), where kind 0=scale, 1=edge, 2=vertex-subset mask.
std::pair<int,Mask> validate_witness(const Graph& g, Mask weights3, Mask tree, int coefficient) {
    if (!is_tree(g,tree) || coefficient<1 || coefficient>K)
        throw std::runtime_error("Invalid witness tree or coefficient");
    const int residual_scale=K-coefficient;
    std::vector<int> residual;
    int total=0;
    std::pair<int,Mask> blocker{-1,0};
    if (coefficient==K) blocker={0,0};
    for (std::size_t i=0; i<g.edges.size(); ++i) {
        int z=2+bool(weights3&(Mask(1)<<i))-coefficient*bool(tree&(Mask(1)<<i));
        if (z<0) throw std::runtime_error("Negative residual coordinate");
        residual.push_back(z); total+=z;
        if (z==0 && (tree&(Mask(1)<<i))) blocker={1,static_cast<Mask>(i)};
    }
    if (total!=residual_scale*(g.n-1)) throw std::runtime_error("Wrong residual total");
    for (Mask vertices=1; vertices<(Mask(1)<<g.n)-1; ++vertices) {
        const int size=pc(vertices);
        if (size<2) continue;
        int inside=0, tree_inside=0;
        for (std::size_t i=0; i<g.edges.size(); ++i) {
            auto [u,v]=g.edges[i];
            if ((vertices&(Mask(1)<<u)) && (vertices&(Mask(1)<<v))) {
                inside+=residual[i];
                tree_inside+=bool(tree&(Mask(1)<<i));
            }
        }
        const int bound=residual_scale*(size-1);
        if (inside>bound) throw std::runtime_error("Residual violates subset inequality");
        if (inside==bound && tree_inside<size-1) blocker={2,vertices};
    }
    if (blocker.first<0) throw std::runtime_error("Witness has no blocking constraint");
    return blocker;
}

struct Counts {
    std::uint64_t graphs=0, candidates=0, feasible=0, trees_considered=0;
    std::map<int,std::uint64_t> coefficients;
};

void write_counts(std::ostream& out, const Counts& c) {
    out << "{\"graphs\":" << c.graphs << ",\"candidates\":" << c.candidates
        << ",\"feasible\":" << c.feasible << ",\"certified_witnesses\":" << c.feasible
        << ",\"trees_considered\":" << c.trees_considered << ",\"coefficient_counts\":{";
    bool first=true;
    for (const auto& [k,v]:c.coefficients) {
        if (!first) out << ',';
        first=false; out << '"' << k << "\":" << v;
    }
    out << "}}";
}

int main() {
    try {
        const auto start=std::chrono::steady_clock::now();
        std::array<std::set<Mask>,8> levels;
        levels[0].insert(0);
        std::array<std::uint64_t,8> extensions{};
        for (int n=1; n<=7; ++n) {
            for (Mask previous:levels[n-1]) {
                Graph old(n-1,previous);
                for (Mask neighbors=0; neighbors<(Mask(1)<<(n-1)); ++neighbors) {
                    ++extensions[n];
                    levels[n].insert(canonical(Graph(n,extend(old,neighbors))));
                }
            }
            std::cout << "unlabeled graphs n=" << n << ": " << levels[n].size() << '\n' << std::flush;
        }
        const auto generated=std::chrono::steady_clock::now();
        std::ofstream catalogue("graph_catalogue.json");
        catalogue << "{\"encoding\":\"bit positions are lexicographic pairs 0<=u<v<n\",\"levels\":[";
        for (int n=0; n<=7; ++n) {
            if (n) catalogue << ',';
            catalogue << "{\"n\":" << n << ",\"masks\":[";
            bool first=true;
            for (Mask mask:levels[n]) { if (!first) catalogue << ','; first=false; catalogue << mask; }
            catalogue << "]}";
        }
        catalogue << "]}\n";
        catalogue.close();

        std::ofstream graphs("graphs.json"), witnesses("witnesses.jsonl");
        graphs << "{\"k\":5,\"vertex_labels\":\"zero based\",\"edge_order\":\"lexicographic\",\"graphs\":[\n";
        std::array<Counts,8> counts;
        std::size_t graph_id=0;
        bool first_graph=true;
        for (int n=5; n<=7; ++n) for (Mask graphmask:levels[n]) {
            Graph g(n,graphmask);
            if (!core_filter(g)) continue;
            const int m=static_cast<int>(g.edges.size());
            const int q=K*(n-1)-2*m;
            assert(q>=0 && q<=m);
            const auto sets=subsets(g);
            const auto all_trees=trees(g,sets);
            Counts local; local.graphs=1;
            combinations(m,q,[&](Mask weights3) {
                ++local.candidates;
                assert(pc(weights3)==q);
                std::vector<int> slacks;
                for (const auto& s:sets) {
                    const int slack=K*s.rank_bound-2*pc(s.edges)-pc(s.edges&weights3);
                    if (slack<0) return;
                    slacks.push_back(slack);
                }
                ++local.feasible;
                bool found=false;
                for (const auto& t:all_trees) {
                    ++local.trees_considered;
                    Rational optimum{K,1};
                    for (int i=0; i<m; ++i) if (t.mask&(Mask(1)<<i))
                        optimum.minimize(2+bool(weights3&(Mask(1)<<i)),1);
                    for (std::size_t i=0; i<sets.size(); ++i)
                        if (t.denominators[i]>0) optimum.minimize(slacks[i],t.denominators[i]);
                    if (optimum.numerator<=0 || optimum.numerator%optimum.denominator) continue;
                    const int coefficient=optimum.numerator/optimum.denominator;
                    const auto blocker=validate_witness(g,weights3,t.mask,coefficient);
                    ++local.coefficients[coefficient];
                    witnesses << "{\"graph_id\":" << graph_id << ",\"weight_three_edge_mask\":" << weights3
                              << ",\"tree_edge_mask\":" << t.mask << ",\"coefficient\":" << coefficient
                              << ",\"blocking_kind\":" << blocker.first << ",\"blocking_index\":" << blocker.second << "}\n";
                    found=true; break;
                }
                if (!found) {
                    std::ostringstream msg;
                    msg << "COUNTEREXAMPLE n=" << n << " graphmask=" << graphmask << " weights3=" << weights3;
                    throw std::runtime_error(msg.str());
                }
            });
            auto& total=counts[n];
            total.graphs+=local.graphs; total.candidates+=local.candidates;
            total.feasible+=local.feasible; total.trees_considered+=local.trees_considered;
            for (auto [k,v]:local.coefficients) total.coefficients[k]+=v;
            if (!first_graph) graphs << ",\n";
            first_graph=false;
            graphs << "{\"id\":" << graph_id << ",\"n\":" << n << ",\"canonical_mask\":" << graphmask
                   << ",\"m\":" << m << ",\"q\":" << q << ",\"edges\":[";
            for (int i=0; i<m; ++i) {
                if (i) graphs << ',';
                graphs << '[' << g.edges[i].first << ',' << g.edges[i].second << ']';
            }
            graphs << "],\"spanning_trees\":" << all_trees.size() << ",\"counts\":";
            write_counts(graphs,local); graphs << '}';
            ++graph_id;
        }
        graphs << "\n]}\n"; graphs.close(); witnesses.close();
        if (!graphs || !witnesses || !catalogue) throw std::runtime_error("Output file write failed");
        const auto finished=std::chrono::steady_clock::now();
        std::ofstream results("results.json");
        results << "{\n\"status\":\"all feasible inputs certified\",\n\"k\":5,\n"
                << "\"generator\":\"all vertex extensions; exact canonicalization over permutations within sorted degree blocks\",\n"
                << "\"arithmetic\":\"integer arithmetic and exact ratio comparisons\",\n"
                << "\"assertions_enabled\":";
#ifdef NDEBUG
        results << "false";
#else
        results << "true";
#endif
        results << ",\n\"unlabeled_graph_counts\":[";
        for (int n=0; n<=7; ++n) { if (n) results << ','; results << levels[n].size(); }
        results << "],\n\"extension_counts\":[";
        for (int n=0; n<=7; ++n) { if (n) results << ','; results << extensions[n]; }
        results << "],\n\"by_vertices\":{";
        Counts total;
        for (int n=5; n<=7; ++n) {
            if (n!=5) results << ',';
            results << '"' << n << "\":"; write_counts(results,counts[n]);
            total.graphs+=counts[n].graphs; total.candidates+=counts[n].candidates;
            total.feasible+=counts[n].feasible; total.trees_considered+=counts[n].trees_considered;
            for (auto [k,v]:counts[n].coefficients) total.coefficients[k]+=v;
            std::cout << "n=" << n << ' '; write_counts(std::cout,counts[n]); std::cout << '\n';
        }
        results << "},\n\"total\":"; write_counts(results,total);
        results << ",\n\"generation_seconds\":" << std::chrono::duration<double>(generated-start).count()
                << ",\n\"audit_seconds\":" << std::chrono::duration<double>(finished-generated).count()
                << ",\n\"total_seconds\":" << std::chrono::duration<double>(finished-start).count() << "\n}\n";
        results.close();
        if (!results) throw std::runtime_error("Results write failed");
        std::cout << "TOTAL "; write_counts(std::cout,total);
        std::cout << "\nruntime_seconds=" << std::chrono::duration<double>(finished-start).count() << '\n';
        return 0;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n'; return 1;
    }
}
