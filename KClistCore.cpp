#include "KClistCore.h"

using namespace std;

#define NLINKS 100000000 //maximum number of edges for memory allocation, will increase if needed

// heap data structure :

struct HashFunc
{
	template<typename T, typename U>
	size_t operator()(const std::pair<T, U>& p) const {
		return std::hash<T>()(p.first) ^ std::hash<U>()(p.second);
	}
};

// Key-value comparison, comparison definition of hash collision, needs to determine whether two custom objects are equal
struct EqualKey {
	template<typename T, typename U>
	bool operator ()(const std::pair<T, U>& p1, const std::pair<T, U>& p2) const {
		return p1.first == p2.first && p1.second == p2.second;
	}
};

typedef struct {
	unsigned key;
	unsigned value;
} keyvalue;

typedef struct {
	unsigned n_max;	// max number of nodes.
	unsigned n;	// number of nodes.
	unsigned* pt;	// pointers to nodes.
	keyvalue* kv; // nodes.
} bheap;

bheap* construct(unsigned n_max) {
	unsigned i;
	bheap* heap = (bheap*)malloc(sizeof(bheap));

	heap->n_max = n_max;
	heap->n = 0;
	heap->pt = (unsigned*)malloc(n_max * sizeof(unsigned));
	for (i = 0; i < n_max; i++) heap->pt[i] = -1;
	heap->kv = (keyvalue*)malloc(n_max * sizeof(keyvalue));
	return heap;
}

inline void swap(bheap* heap, unsigned i, unsigned j) {
	keyvalue kv_tmp = heap->kv[i];
	unsigned pt_tmp = heap->pt[kv_tmp.key];
	heap->pt[heap->kv[i].key] = heap->pt[heap->kv[j].key];
	heap->kv[i] = heap->kv[j];
	heap->pt[heap->kv[j].key] = pt_tmp;
	heap->kv[j] = kv_tmp;
}

inline void bubble_up(bheap* heap, unsigned i) {
	unsigned j = (i - 1) / 2;
	while (i > 0) {
		if (heap->kv[j].value > heap->kv[i].value) {
			swap(heap, i, j);
			i = j;
			j = (i - 1) / 2;
		}
		else break;
	}
}

inline void bubble_down(bheap* heap) {
	unsigned i = 0, j1 = 1, j2 = 2, j;
	while (j1 < heap->n) {
		j = ((j2 < heap->n) && (heap->kv[j2].value < heap->kv[j1].value)) ? j2 : j1;
		if (heap->kv[j].value < heap->kv[i].value) {
			swap(heap, i, j);
			i = j;
			j1 = 2 * i + 1;
			j2 = j1 + 1;
			continue;
		}
		break;
	}
}

inline void insert(bheap* heap, keyvalue kv) {
	heap->pt[kv.key] = (heap->n)++;
	heap->kv[heap->n - 1] = kv;
	bubble_up(heap, heap->n - 1);
}

inline void update(bheap* heap, unsigned key) {
	unsigned i = heap->pt[key];
	if (i != -1) {
		((heap->kv[i]).value)--;
		bubble_up(heap, i);
	}
}

inline keyvalue popmin(bheap* heap) {
	keyvalue min = heap->kv[0];
	heap->pt[min.key] = -1;
	heap->kv[0] = heap->kv[--(heap->n)];
	heap->pt[heap->kv[0].key] = 0;
	bubble_down(heap);
	return min;
}

// graph datastructure:

typedef struct {
	unsigned s;
	unsigned t;
} edge;

typedef struct {
	//edge list structure:
	unsigned n; //number of nodes
	unsigned e; //number of edges
	unsigned n2; //number of nodes with core value larger than one
	unsigned e2; //number of edges between nodes with core value larger than one
	edge* edges;//list of edges

	//to compute a degeneracy ordering:
	unsigned* d0; //degrees
	unsigned* cd0; //cumulative degree: (start with 0) length=dim+1
	unsigned* adj0; //list of neighbors
	unsigned* rank; //degeneracy rankings of nodes
	unsigned* map;//map[newlabel]=oldlabel
	unsigned core; //core number of the graph

	//truncated neighborhoods:
	unsigned* d; //truncated degrees
	unsigned* cd; //cumulative degree: (start with 0) length=dim+1
	unsigned* adj; //list of neighbors with higher rank

	unsigned kmax;
} sparse;

unordered_map<pair<int, int>, int, HashFunc, EqualKey>mpall;

//compute the maximum of three unsigned
inline unsigned max3(unsigned a, unsigned b, unsigned c) {
	a = (a > b) ? a : b;
	return (a > c) ? a : c;
}

//reading the edgelist from file
sparse* readedgelist(char* edgelist) {
	unsigned e1 = NLINKS;
	sparse* g = (sparse*)malloc(sizeof(sparse));
	FILE* file;

	g->n = 0;
	g->e = 0;
	file = fopen(edgelist, "r");
	g->edges = (edge*)malloc(e1 * sizeof(edge));
	while (fscanf(file, "%u %u", &(g->edges[g->e].s), &(g->edges[g->e].t)) == 2) {
		g->n = max3(g->n, g->edges[g->e].s, g->edges[g->e].t);
		if (g->e++ == e1) {
			e1 += NLINKS;
			g->edges = (edge*)realloc(g->edges, e1 * sizeof(edge));  //g->edges 存取边 [].s/t 编号对应s和t
		}
	}
	fclose(file);
	g->n++;

	g->edges = (edge*)realloc(g->edges, g->e * sizeof(edge));

	return g;
}

sparse* getEdgeList(std::vector<std::vector<int>> adj, vector<bool>selected) {
	unsigned e1 = NLINKS;
	sparse* g = (sparse*)malloc(sizeof(sparse));


	g->n = 0;
	g->e = 0;

	g->edges = (edge*)malloc(e1 * sizeof(edge));

	unordered_map<pair<int, int>, int, HashFunc, EqualKey>mp;

	for (int i = 0; i < adj.size(); i++)
	{
		if (selected[i] == false)continue;
		for (int j = 0; j < adj[i].size(); j++)
		{
			int s = i, t = adj[i][j];
			if (selected[t] == false)continue;
			if (mp.find({ s,t }) == mp.end())
			{
				mp.insert({ { t,s }, 1 });
				g->edges[g->e].s = s, g->edges[g->e].t = t;
				g->n = max3(g->n, g->edges[g->e].s, g->edges[g->e].t);
				if (g->e++ == e1) {
					e1 += NLINKS;
					g->edges = (edge*)realloc(g->edges, e1 * sizeof(edge)); 
				}
			}
		}
	}
	g->n++;

	g->edges = (edge*)realloc(g->edges, g->e * sizeof(edge));

	return g;
}

sparse* getEdgeList_hash(std::vector<std::vector<int>> adj, vector<bool>selected) {
	unsigned e1 = NLINKS;
	sparse* g = (sparse*)malloc(sizeof(sparse));
	g->n = 0;
	g->e = 0;
	g->edges = (edge*)malloc(e1 * sizeof(edge));

	unordered_map<pair<int, int>, int, HashFunc, EqualKey>mp;

	for (int i = 0; i < adj.size(); i++)
	{
		if (selected[i] == false)continue;
		for (int j = 0; j < adj[i].size(); j++)
		{
			int s = i, t = adj[i][j];
			if (selected[t] == false)continue;
			if (mp.find({ s,t }) == mp.end())
			{
				mp.insert({ { t,s }, 1 });
				if (mpall.find({ s,t }) == mpall.end())
				{
					mpall.insert({ { s,t }, 1 });
					mpall.insert({ { t,s }, 1 });
				}		
				g->edges[g->e].s = s, g->edges[g->e].t = t;
				g->n = max3(g->n, g->edges[g->e].s, g->edges[g->e].t);
				if (g->e++ == e1) {
					e1 += NLINKS;
					g->edges = (edge*)realloc(g->edges, e1 * sizeof(edge));  
				}
			}
		}
	}
	g->n++;

	g->edges = (edge*)realloc(g->edges, g->e * sizeof(edge));

	return g;
}

//Building the graph structure
void mkgraph(sparse* g) {
	unsigned i;
	g->d0 = (unsigned*)calloc(g->n, sizeof(unsigned));

	for (i = 0; i < g->e; i++) {
		g->d0[g->edges[i].s]++;
		g->d0[g->edges[i].t]++;
	}   // degree
	g->cd0 = (unsigned*)malloc((g->n + 1) * sizeof(unsigned)); // sum degree of d[0] - d[i-1]
	g->cd0[0] = 0;
	g->kmax = 0;
	for (i = 1; i < g->n + 1; i++) {
		g->kmax = max(g->kmax, g->d0[i - 1]);
		g->cd0[i] = g->cd0[i - 1] + g->d0[i - 1];
		g->d0[i - 1] = 0;
	}
	/*
	for (int i = 0; i < g->n + 1; i++)
	{
		printf("%d:%d\n", i, g->cd0[i]);

		printf("%d:%d\n", i, g->d0[i]);
	}*/

	g->adj0 = (unsigned*)malloc(2 * g->e * sizeof(unsigned));

	for (i = 0; i < g->e; i++) {
		g->adj0[g->cd0[g->edges[i].s] + g->d0[g->edges[i].s]++] = g->edges[i].t;
		g->adj0[g->cd0[g->edges[i].t] + g->d0[g->edges[i].t]++] = g->edges[i].s;
	}
	/*
	for (int i = 0; i < g->e; i++)
	{
		printf("%d:%d\n", i, g->adj0[i]);
	}*/

}

//Building the heap structure with (key,value)=(node,degree) for each node
bheap* mkheap(sparse* g) {
	unsigned i;
	keyvalue kv;
	bheap* heap = construct(g->n);
	for (i = 0; i < g->n; i++) {
		kv.key = i;
		kv.value = g->d0[i];
		insert(heap, kv);
	}
	return heap;
}

void freeheap(bheap* heap) {
	free(heap->pt);
	free(heap->kv);
	free(heap);
}

void getrankloop(sparse* g, unsigned kmax)
{
	g->rank = (unsigned*)malloc(g->n * sizeof(unsigned));  
	g->rank = (unsigned*)malloc(g->n * sizeof(unsigned)); 
	for (int i = 0; i < g->n; i++)
	{
		g->rank[i] = i;
	}
	g->n2 = g->n;
}

//computing degeneracy ordering and core value
void kcore(sparse* g, unsigned kmax) {
	unsigned i, j, r = 0, n = g->n, k = kmax - 1;  
	keyvalue kv;   
	unsigned c = 0;//the core number
	bheap* heap = mkheap(g);
	g->rank = (unsigned*)malloc(g->n * sizeof(unsigned));  
	g->map = (unsigned*)malloc(g->n * sizeof(unsigned));
	for (i = 0; i < g->n; i++) {
		kv = popmin(heap);
		if (kv.value > c) {
			c = kv.value;
		}
		if (c < k) {//remove node with core value less than kmax-1
			g->rank[kv.key] = -1;
			n--;
		}
		else {
			g->map[n - (++r)] = kv.key;
			g->rank[kv.key] = n - r;
		}
		for (j = g->cd0[kv.key]; j < g->cd0[kv.key + 1]; j++) {
			update(heap, g->adj0[j]);
		}
	}
	/*
	for (int i = 0; i < n; i++)
	{
		printf("%d:%d\n", i, g->rank[i]);
	}*/

	freeheap(heap);
	/*free(g->d0);
	free(g->cd0);
	free(g->adj0);*/
	g->core = c;
	g->n2 = n;
}

void relabelnodes(sparse* g) {
	unsigned i, j, source, target;
	j = 0;
	for (i = 0; i < g->e; i++) {
		source = g->rank[g->edges[i].s];
		target = g->rank[g->edges[i].t];
		if (source == -1 || target == -1) {
			continue;
		}
		if (source < target) {
			g->edges[j].s = target;
			g->edges[j++].t = source;
		}
		else {
			g->edges[j].s = source;
			g->edges[j++].t = target;
		}
	}
	g->e2 = j;
	g->edges = (edge*)realloc(g->edges, g->e2 * sizeof(edge));
	/*
	for (int i = 0; i < g->e2; i++)
	{
		printf("%d %d\n", g->edges[i].s, g->edges[i].t);
	}*/
}

//for future use in qsort
int cmpfunc(const void* a, const void* b) {
	if (*(unsigned*)a > *(unsigned*)b) {
		return 1;
	}
	return -1;
}


//Building the special graph structure
void mkspecial(sparse* g) {
	unsigned i;
	g->d = (unsigned*)calloc(g->n2, sizeof(unsigned));

	//g->d 
	for (i = 0; i < g->e2; i++) {
		g->d[g->edges[i].s]++;
	}



	g->cd = (unsigned*)malloc((g->n2 + 1) * sizeof(unsigned));
	g->cd[0] = 0;
	for (i = 1; i < g->n2 + 1; i++) {
		g->cd[i] = g->cd[i - 1] + g->d[i - 1];
		g->d[i - 1] = 0;
	}
	/*
	for (int i = 0; i < g->n2 + 1; i++)
	{
		printf("%d:%d %d\n", i, g->cd[i], g->d[i]);
	}*/

	g->adj = (unsigned*)malloc((g->e2) * sizeof(unsigned));

	for (i = 0; i < g->e2; i++) {
		g->adj[g->cd[g->edges[i].s] + g->d[g->edges[i].s]++] = g->edges[i].t;
	}
	/*
	for (int i = 0; i < g->e2; i++)
	{
		printf("%d:%d\n", i, g->adj[i]);
	}*/

	for (i = 0; i < g->n2; i++) {
		qsort(&g->adj[g->cd[i]], g->d[i], sizeof(unsigned), cmpfunc);
	}
	//free(g->edges); Can be freed if node parallelisation is used instead of edge
}




void freesparse(sparse* g) {
	free(g->edges);
	free(g->map);
	free(g->rank);
	free(g->d);
	free(g->cd);
	free(g->adj);
	free(g->d0);
	free(g->cd0);
	free(g->adj0);
	free(g);
}

//store the intersection of list1 and list2 in list3 and return the size of list3 (the 3 lists are sorted)
inline unsigned merging(unsigned* list1, unsigned s1, unsigned* list2, unsigned s2, unsigned* list3) {
	unsigned i = 0, j = 0, s3 = 0;
	unsigned x = list1[0], y = list2[0];
	while (i < s1 && j < s2) {
		if (x < y) {
			x = list1[++i];
			continue;
		}
		if (y < x) {
			y = list2[++j];
			continue;
		}
		list3[s3++] = x;
		x = list1[++i];
		y = list2[++j];
	}
	return s3;
}


//store the intersection of list1 and list2 in list3 and return the size of list3 (the 3 lists are sorted)
inline unsigned merging(unsigned* list1, unsigned s1, unsigned* list2, unsigned s2, unsigned* list3, unsigned w) {
	unsigned i = 0, j = 0, s3 = 0;
	unsigned x = list1[0], y = list2[0];
	while (i < s1 && j < s2) {
		if (x < y) {
			x = list1[++i];
			continue;
		}
		if (y < x) {
			y = list2[++j];
			continue;
		}
		if (x != w)list3[s3++] = x;
		x = list1[++i];
		y = list2[++j];
	}
	return s3;
}

void generateCombinations(unsigned* merge, int mergeSize, sparse* g, unsigned* ck, int ckSize, std::vector<std::vector<int>>& out, vector<int>temp, int start, int k) {
	if (k == 0) {
		out.push_back(temp);
		return;
	}

	for (int i = start; i < mergeSize; ++i) {
		temp.push_back(g->map[(int)merge[i]]);
		generateCombinations(merge, mergeSize, g, ck, ckSize, out, temp, i + 1, k - 1);
		temp.pop_back();
	}
}

//the recursion to compute all possible intersections
void recursion(unsigned kmax, unsigned k, unsigned* merge, unsigned* size, sparse* g, unsigned* ck, unsigned long long* nck, std::vector<std::vector<int>>& clique_out, vector<vector<int>>& v2cliques) {
	unsigned t = (k - 3) * g->core, t2 = t + g->core;
	unsigned i, u;

	if (size[k - 3] < kmax - k) {//stop if we already know k-cliques cannot be formed
		return;
	}

	if (k == kmax) {//increasing the k-clique degrees

		if (size[k - 3] <= 0)return;
		int id = clique_out.size();
		vector<int>temp;
		for (int i = 0; i < kmax - 1; i++)
		{
			temp.push_back(g->map[ck[i]]);
			for (int j = 0; j < size[k - 3]; j++)
			{
				v2cliques[g->map[ck[i]]].push_back(id + j);
			}
		}
		for (int i = 0; i < size[k - 3]; i++)
		{
			temp.push_back(g->map[merge[t + i]]);
			v2cliques[g->map[merge[t + i]]].push_back(id + i);
			clique_out.push_back(temp);
			temp.pop_back();
		}
		/*
		for (int j = 1; j <= size[k - 3]; ++j) {
			generateCombinations(merge + t, size[k - 3], g, ck, kmax - 1, clique_out, temp, 0, j);
		}*/
		return;
	}

	for (i = 0; i < size[k - 3]; i++) {
		ck[k - 1] = merge[t + i];
		//printf("We add %d(ck=%d)\n", g->map[merge[t + i]],ck[k-1]);
		size[k - 2] = merging(&g->adj[g->cd[ck[k - 1]]], g->d[ck[k - 1]], &merge[t], size[k - 3], &merge[t2]);
		recursion(kmax, k + 1, merge, size, g, ck, nck, clique_out, v2cliques);
	}
}

//one pass over all k-cliques
unsigned long long* onepass(sparse* g, unsigned kmax, std::vector<std::vector<int>>& clique_out, vector<vector<int>>& v2cliques) {
	unsigned e, i;
	unsigned* merge, * size, * ck;
	unsigned long long* nck;

	merge = (unsigned*)malloc((kmax - 2) * g->core * sizeof(unsigned));
	size = (unsigned*)malloc((kmax - 2) * sizeof(unsigned));
	ck = (unsigned*)malloc(kmax * sizeof(unsigned));
	nck = (unsigned long long*)calloc(g->n2, sizeof(unsigned long long));

	for (e = 0; e < g->e2; e++) {
		ck[0] = g->edges[e].s;
		ck[1] = g->edges[e].t;
		size[0] = merging(&(g->adj[g->cd[ck[0]]]), g->d[ck[0]], &(g->adj[g->cd[ck[1]]]), g->d[ck[1]], merge);
		recursion(kmax, 3, merge, size, g, ck, nck, clique_out, v2cliques);
	}
	free(merge);
	free(size);
	free(ck);

	return nck;
}


void onepassforloop(sparse* g, unsigned kmax, std::vector<std::vector<int>>& clique_out, vector<vector<int>>& v2cliques) {
	unsigned* merge, * size, * ck;

	merge = (unsigned*)malloc(g->kmax * sizeof(unsigned));
	size = (unsigned*)malloc(g->kmax * sizeof(unsigned));
	ck = (unsigned*)malloc(4 * sizeof(unsigned));

	for (int i = 0; i < g->n2; i++)
	{
		vector<int>temp;
		int degi = g->d[i];
		ck[0] = i;
		temp.push_back(ck[0]);
		for (int u = 0; u < g->d[i]; u++)
		{
			ck[1] = g->adj[g->cd[i] + u];
			temp.push_back(ck[1]);
			for (int v = u + 1; v < g->d[i]; v++)
			{
				ck[2] = g->adj[g->cd[i] + v];
				temp.push_back(ck[2]);
				size[0] = merging(&(g->adj0[g->cd0[ck[1]]]), g->d0[ck[1]], &(g->adj0[g->cd0[ck[2]]]), g->d0[ck[2]], merge, ck[0]);
				for (int w = 0; w < size[0]; w++)
				{
					if (g->rank[merge[w]] > g->rank[ck[0]])continue;
					ck[3] = merge[w];
					temp.push_back(ck[3]);
					clique_out.push_back(temp);
					int id = clique_out.size() - 1;
					for (int zz = 0; zz <= 3; zz++)
					{
						v2cliques[ck[zz]].push_back(id);
					}
					temp.pop_back();
				}
				temp.pop_back();
				//for(int )
			}
			temp.pop_back();
		}
	}
	free(merge);
	free(size);
	free(ck);


	return;
}

void onepassfordoubletri(sparse* g, unsigned kmax, std::vector<std::vector<int>>& clique_out, vector<vector<int>>& v2cliques)
{
	unsigned* merge, * size, * ck;

	merge = (unsigned*)malloc(g->kmax * sizeof(unsigned));
	size = (unsigned*)malloc(g->kmax * sizeof(unsigned));
	ck = (unsigned*)malloc(4 * sizeof(unsigned));


	for (int e = 0; e < g->e2; e++) {
		vector<int>temp;
		ck[0] = g->edges[e].s;
		ck[1] = g->edges[e].t;
		temp.push_back(ck[0]);
		temp.push_back(ck[1]);
		size[0] = merging(&(g->adj0[g->cd0[ck[0]]]), g->d0[ck[0]], &(g->adj0[g->cd0[ck[1]]]), g->d0[ck[1]], merge);

		for (int w = 0; w < size[0]; w++)
		{
			ck[2] = merge[w];
			temp.push_back(ck[2]);
			for (int ww = w + 1; ww < size[0]; ww++)
			{
				ck[3] = merge[ww];
				temp.push_back(ck[3]);
				clique_out.push_back(temp);
				int id = clique_out.size() - 1;
				for (int zz = 0; zz <= 3; zz++)
				{
					v2cliques[ck[zz]].push_back(id);
				}
				temp.pop_back();
			}
			temp.pop_back();
		}

	}
	free(merge);
	free(size);
	free(ck);


	return;
}

void dfsstar(unsigned* merge, unsigned* size, int start, int k, vector<int>& combination, std::vector<std::vector<int>>& clique_out, vector<vector<int>>& v2cliques) {
	if (k == 0) {
		clique_out.push_back(combination);
		int id = clique_out.size() - 1;
		for (auto t : combination)
		{
			v2cliques[t].push_back(id);
		}
		return;
	}

	for (int i = start; i < size[0]; ++i) {
		combination.push_back(merge[i]);
		dfsstar(merge, size, i + 1, k - 1, combination, clique_out, v2cliques);
		combination.pop_back();
	}
}

void dfsctri(unsigned* merge, unsigned* size, int start, int k, vector<int>& combination, std::vector<std::vector<int>>& clique_out, vector<vector<int>>& v2cliques) {
	if (k == 0) {
		int tx = combination[1], ty = combination[2], tz = combination[3];
		if (mpall[{tx,ty}])
		{
			clique_out.push_back(combination);
			int id = clique_out.size() - 1;
			for (auto t : combination)
			{
				v2cliques[t].push_back(id);
			}
		}
		if (mpall.find({ tx,tz }) != mpall.end())
		{
			clique_out.push_back(combination);
			int id = clique_out.size() - 1;
			for (auto t : combination)
			{
				v2cliques[t].push_back(id);
			}
		}
		if (mpall.find({ ty,tz }) != mpall.end())
		{
			clique_out.push_back(combination);
			int id = clique_out.size() - 1;
			for (auto t : combination)
			{
				v2cliques[t].push_back(id);
			}
		}
		return;
	}

	for (int i = start; i < size[0]; ++i) {
		combination.push_back(merge[i]);
		dfsctri(merge, size, i + 1, k - 1, combination, clique_out, v2cliques);
		combination.pop_back();
	}
}

void onepassforstar(sparse* g, unsigned kmax, std::vector<std::vector<int>>& clique_out, vector<vector<int>>& v2cliques)
{
	unsigned* merge, * size, * ck;

	merge = (unsigned*)malloc(g->n2 * sizeof(unsigned));
	size = (unsigned*)malloc(g->n2 * sizeof(unsigned));
	ck = (unsigned*)malloc(g->n2 * sizeof(unsigned));


	for (int i = 0; i < g->n2; i++)
	{
		vector<int>temp;
		ck[0] = i;
		temp.push_back(ck[0]);
		size[0] = 0;
		for (int u = 0; u < g->d0[i]; u++)
		{
			merge[size[0]++] = g->adj0[g->cd0[i] + u];
		}
		dfsstar(merge, size, 0, kmax, temp, clique_out, v2cliques);
	}
	free(merge);
	free(size);
	free(ck);

	return;
}


void onepassforctri(sparse* g, unsigned kmax, std::vector<std::vector<int>>& clique_out, vector<vector<int>>& v2cliques)
{
	unsigned* merge, * size, * ck;

	merge = (unsigned*)malloc(g->n2 * sizeof(unsigned));
	size = (unsigned*)malloc(g->n2 * sizeof(unsigned));
	ck = (unsigned*)malloc(g->n2 * sizeof(unsigned));


	for (int i = 0; i < g->n2; i++)
	{
		vector<int>temp;
		ck[0] = i;
		temp.push_back(ck[0]);
		size[0] = 0;
		for (int u = 0; u < g->d0[i]; u++)
		{
			merge[size[0]++] = g->adj0[g->cd0[i] + u];
		}
		dfsctri(merge, size, 0, 3, temp, clique_out, v2cliques);
	}
	free(merge);
	free(size);
	free(ck);

	return;
}

void onepassfortritri(sparse* g, unsigned kmax, std::vector<std::vector<int>>& clique_out, vector<vector<int>>& v2cliques)
{
	unsigned* merge, * size, * ck;

	merge = (unsigned*)malloc(g->kmax * sizeof(unsigned));
	size = (unsigned*)malloc(g->kmax * sizeof(unsigned));
	ck = (unsigned*)malloc(5 * sizeof(unsigned));


	for (int e = 0; e < g->e2; e++) {
		vector<int>temp;
		ck[0] = g->edges[e].s;
		ck[1] = g->edges[e].t;
		temp.push_back(ck[0]);
		temp.push_back(ck[1]);
		size[0] = merging(&(g->adj0[g->cd0[ck[0]]]), g->d0[ck[0]], &(g->adj0[g->cd0[ck[1]]]), g->d0[ck[1]], merge);

		for (int w = 0; w < size[0]; w++)
		{
			ck[2] = merge[w];
			temp.push_back(ck[2]);
			for (int ww = w + 1; ww < size[0]; ww++)
			{
				ck[3] = merge[ww];
				temp.push_back(ck[3]);
				for (int www = ww + 1; www < size[0]; www++)
				{
					ck[4] = merge[www];
					temp.push_back(ck[4]);
					clique_out.push_back(temp);
					int id = clique_out.size() - 1;
					for (int zz = 0; zz <= 4; zz++)
					{
						v2cliques[ck[zz]].push_back(id);
					}
					temp.pop_back();
				}
				temp.pop_back();
			}
			temp.pop_back();
		}

	}
	free(merge);
	free(size);
	free(ck);


	return;
}


//printing the k-clique degree of each node
unsigned long long printckdeg(unsigned long long* nck, sparse* g, char* ckdeg) {
	unsigned i;
	unsigned long long tot = 0;
	FILE* file = fopen(ckdeg, "w");
	for (i = 0; i < g->n2; i++) {
		fprintf(file, "%u %llu\n", g->map[i], nck[i]);
		tot += nck[i];
	}
	fclose(file);
	return tot;
}

// heap data structure :

typedef struct {
	unsigned key;
	unsigned long long value;
} keyvalueLLU;

typedef struct {
	unsigned n_max;// max number of nodes.
	unsigned n;// number of nodes.
	unsigned* pt;// pointers to nodes.
	keyvalueLLU* kv;// (node,nck)
} bheapLLU;

bheapLLU* constructLLU(unsigned n_max) {
	unsigned i;
	bheapLLU* heap = (bheapLLU*)malloc(sizeof(bheapLLU));

	heap->n_max = n_max;
	heap->n = 0;
	heap->pt = (unsigned*)malloc(n_max * sizeof(unsigned));
	for (i = 0; i < n_max; i++) heap->pt[i] = -1;
	heap->kv = (keyvalueLLU*)malloc(n_max * sizeof(keyvalueLLU));
	return heap;
}

inline void swapLLU(bheapLLU* heap, unsigned i, unsigned j) {
	keyvalueLLU kv_tmp = heap->kv[i];
	unsigned pt_tmp = heap->pt[kv_tmp.key];
	heap->pt[heap->kv[i].key] = heap->pt[heap->kv[j].key];
	heap->kv[i] = heap->kv[j];
	heap->pt[heap->kv[j].key] = pt_tmp;
	heap->kv[j] = kv_tmp;
}

inline void bubble_upLLU(bheapLLU* heap, unsigned i) {
	unsigned j = (i - 1) / 2;
	while (i > 0) {
		if (heap->kv[j].value > heap->kv[i].value) {
			swapLLU(heap, i, j);
			i = j;
			j = (i - 1) / 2;
		}
		else break;
	}
}

inline void bubble_downLLU(bheapLLU* heap) {
	unsigned i = 0, j1 = 1, j2 = 2, j;
	while (j1 < heap->n) {
		j = ((j2 < heap->n) && (heap->kv[j2].value < heap->kv[j1].value)) ? j2 : j1;
		if (heap->kv[j].value < heap->kv[i].value) {
			swapLLU(heap, i, j);
			i = j;
			j1 = 2 * i + 1;
			j2 = j1 + 1;
			continue;
		}
		break;
	}
}

inline void insertLLU(bheapLLU* heap, keyvalueLLU kv) {
	heap->pt[kv.key] = (heap->n)++;
	heap->kv[heap->n - 1] = kv;
	bubble_upLLU(heap, heap->n - 1);
}

inline void updateLLU(bheapLLU* heap, unsigned key, unsigned long long delta) {
	unsigned i = heap->pt[key];
	if (i != -1) {
		((heap->kv[i]).value) -= delta;
		bubble_upLLU(heap, i);
	}
}

inline keyvalueLLU popminLLU(bheapLLU* heap) {
	keyvalueLLU min = heap->kv[0];
	heap->pt[min.key] = -1;
	heap->kv[0] = heap->kv[--(heap->n)];
	heap->pt[heap->kv[0].key] = 0;
	bubble_downLLU(heap);
	return min;
}

//Building the heap structure with (key,value)=(node,k-clique degree) for each node
bheapLLU* mkheapLLU(unsigned long long* nck, unsigned n) {
	unsigned i;
	keyvalueLLU kv;
	bheapLLU* heap = constructLLU(n);
	for (i = 0; i < n; i++) {
		kv.key = i;
		kv.value = nck[i];
		insertLLU(heap, kv);
	}
	return heap;
}

void freeheapLLU(bheapLLU* heap) {
	free(heap->pt);
	free(heap->kv);
	free(heap);
}

typedef struct {
	unsigned n;//number of nodes;
	keyvalueLLU* ic;// (node,k-clique core) in k-clique core ordering
	unsigned long long ckcore;// k-clique core number of the graph
	unsigned size;// size of the k-clique dense subgraph
	double rho;// edge density of the denses subgraph found
	double ckrho;// k-clique density
} densest;

//Building the graph structure
void mkgraph2(sparse* g) {
	unsigned i;
	g->d0 = (unsigned*)calloc(g->n2, sizeof(unsigned));

	for (i = 0; i < g->e2; i++) {
		g->d0[g->edges[i].s]++;
		g->d0[g->edges[i].t]++;
	}
	g->cd0 = (unsigned*)malloc((g->n2 + 1) * sizeof(unsigned));
	g->cd0[0] = 0;
	for (i = 1; i < g->n2 + 1; i++) {
		g->cd0[i] = g->cd0[i - 1] + g->d0[i - 1];
		g->d0[i - 1] = 0;
	}

	g->adj0 = (unsigned*)malloc(2 * g->e2 * sizeof(unsigned));

	for (i = 0; i < g->e2; i++) {
		g->adj0[g->cd0[g->edges[i].s] + g->d0[g->edges[i].s]++] = g->edges[i].t;
		g->adj0[g->cd0[g->edges[i].t] + g->d0[g->edges[i].t]++] = g->edges[i].s;
	}

	for (i = 0; i < g->n2; i++) {
		qsort(&g->adj0[g->cd0[i]], g->d0[i], sizeof(unsigned), cmpfunc);
	}

}


//one pass over all k-cliques
void oneshortpass(sparse* g, bheapLLU* heap, unsigned kmax, unsigned u, unsigned long long* nck) {
	unsigned i, v, s = 0;
	static unsigned* merge = NULL, * size, * ck, * adj;
	if (merge == NULL) {
		merge = (unsigned*)malloc((kmax - 2) * g->core * sizeof(unsigned));
		size = (unsigned*)malloc((kmax - 2) * sizeof(unsigned));
		ck = (unsigned*)malloc(kmax * sizeof(unsigned));
		adj = (unsigned*)malloc(g->n2 * sizeof(unsigned));
	}

	for (i = g->cd0[u]; i < g->cd0[u + 1]; i++) {
		v = g->adj0[i];
		if (heap->pt[v] != -1) {
			adj[s++] = v;
		}
	}
	std::vector<std::vector<int>>clique_out;
	vector<vector<int>>v2cliques;

	if (s > kmax - 2) {
		ck[0] = u;
		for (i = 0; i < s; i++) {
			ck[1] = adj[i];
			size[0] = merging(adj, s, &(g->adj[g->cd[ck[1]]]), g->d[ck[1]], merge);
			recursion(kmax, 3, merge, size, g, ck, nck, clique_out, v2cliques);
		}
	}

}

densest* kcliquecore(unsigned kmax, unsigned long long* nck, unsigned long long ncktot, sparse* g) {
	unsigned i, j, u, v;
	keyvalueLLU kv;
	unsigned long long c = 0;//the k-clique core number
	bheapLLU* heap = mkheapLLU(nck, g->n2);
	unsigned e = g->e2; //number of edges left
	densest* ds = (densest*)malloc(sizeof(densest));
	ds->ic = (keyvalueLLU*)calloc(g->n2, sizeof(keyvalueLLU));
	ds->n = g->n2;
	double ckrho_tmp;
	unsigned size_tmp;
	unsigned long long* nck2 = (unsigned long long*)calloc(g->n2, sizeof(unsigned long long));
	unsigned long long ncktot2;

	ds->ckrho = ((double)ncktot) / ((double)g->n2);
	ds->size = g->n2;
	ds->rho = ((double)(2 * e)) / ((double)(g->n2 * (g->n2 - 1)));
	for (i = 0; i < g->n2 - 1; i++) {
		kv = popminLLU(heap);
		u = kv.key;
		if (kv.value > c) {
			c = kv.value;
		}
		ds->ic[i].key = u;
		ds->ic[i].value = c;
		oneshortpass(g, heap, kmax, u, nck2);
		ncktot2 = 0;
		for (j = g->cd0[u]; j < g->cd0[u + 1]; j++) {
			v = g->adj0[j];
			if (heap->pt[v] != -1) {
				e--;
				updateLLU(heap, v, nck2[v]);
				ncktot2 += nck2[v];
				nck2[v] = 0;
			}
		}
		ncktot -= ncktot2 / ((unsigned long long)(kmax - 1));//kcliques are counted kmax-1 times
		size_tmp = g->n2 - i - 1;
		ckrho_tmp = ((double)ncktot) / ((double)size_tmp);
		if (ckrho_tmp > ds->ckrho) {
			ds->ckrho = ckrho_tmp;
			ds->size = size_tmp;
			ds->rho = ((double)(2 * e)) / ((double)((size_tmp) * (size_tmp - 1)));
		}
	}
	kv = popminLLU(heap);
	ds->ic[i].key = kv.key;
	ds->ic[i].value = c;
	ds->ckcore = c;
	freeheapLLU(heap);
	return ds;
}

void printdensest(densest* ds, sparse* g, char* ckcore, char* ckdens) {
	FILE* file;
	unsigned i;

	file = fopen(ckcore, "w");
	for (i = 0; i < ds->n; i++) {
		fprintf(file, "%u %llu\n", g->map[ds->ic[i].key], ds->ic[i].value);
	}
	fclose(file);

	file = fopen(ckdens, "w");
	fprintf(file, "%lf %lf %u", ds->ckrho, ds->rho, ds->size);
	for (i = ds->n - 1; i > ds->n - ds->size - 1; i--) {
		fprintf(file, " %u", g->map[ds->ic[i].key]);
	}
	fprintf(file, "\n");
	fclose(file);
}

std::vector<std::vector<int>>get_kloops(std::vector<std::vector<int>> adj, vector<vector<int>>& v2cliques, unsigned k_max, vector<bool>selected)
{
	sparse* g;
	unsigned kmax = k_max;
	unsigned long long* nck;
	unsigned long long tot;
	densest* ds;
	g = getEdgeList(adj, selected);

	mkgraph(g);
	//kcore(g, kmax);
	getrankloop(g, kmax);
	relabelnodes(g);
	mkspecial(g);

	std::vector<std::vector<int>>clique_out;
	onepassforloop(g, kmax, clique_out, v2cliques);

	free(g->edges);
	free(g->rank);
	free(g->d);
	free(g->cd);
	free(g->adj);
	free(g->d0);
	free(g->cd0);
	free(g->adj0);
	free(g);
	return clique_out;
}

std::vector<std::vector<int>>get_kcliques(std::vector<std::vector<int>> adj, vector<vector<int>>& v2cliques, unsigned k_max, vector<bool>selected)
{
	sparse* g;
	unsigned kmax = k_max;
	unsigned long long* nck;
	unsigned long long tot;
	densest* ds;
	g = getEdgeList(adj, selected);

	mkgraph(g);
	kcore(g, kmax);
	relabelnodes(g);
	mkspecial(g);

	std::vector<std::vector<int>>clique_out;
	nck = onepass(g, kmax, clique_out, v2cliques);
	/*for (int i = 0; i < clique_out.size(); i++)
	{
		printf("%d-kcliques is found:",i);
		for (int j = 0; j < clique_out[i].size(); j++)
		{
			printf("%d ", clique_out[i][j]);
		}
		printf("\n");
	}*/

	freesparse(g);
	return clique_out;
}

std::vector<std::vector<int>>get_kstar(std::vector<std::vector<int>> adj, vector<vector<int>>& v2cliques, unsigned k_max, vector<bool>selected)
{
	sparse* g;
	unsigned kmax = k_max;
	unsigned long long* nck;
	unsigned long long tot;
	g = getEdgeList(adj, selected);

	mkgraph(g);
	getrankloop(g, kmax);
	relabelnodes(g);
	mkspecial(g);

	std::vector<std::vector<int>>clique_out;
	onepassforstar(g, kmax, clique_out, v2cliques);
	freesparse(g);
	return clique_out;
}

std::vector<std::vector<int>>get_ctriangle(std::vector<std::vector<int>> adj, std::vector<std::vector<int>>& v2cliques, unsigned k_max, std::vector<bool>selected)
{
	sparse* g;
	unsigned kmax = k_max;
	unsigned long long* nck;
	unsigned long long tot;
	g = getEdgeList_hash(adj, selected);

	mkgraph(g);
	getrankloop(g, kmax);
	relabelnodes(g);
	mkspecial(g);

	std::vector<std::vector<int>>clique_out;
	onepassforctri(g, kmax, clique_out, v2cliques);
	free(g->edges);
	free(g->rank);
	free(g->d);
	free(g->cd);
	free(g->adj);
	free(g->d0);
	free(g->cd0);
	free(g->adj0);
	free(g);
	return clique_out;
}

std::vector<std::vector<int>>get_doubletriangle(std::vector<std::vector<int>> adj, std::vector<std::vector<int>>& v2cliques, unsigned k_max, std::vector<bool>selected)
{
	sparse* g;
	unsigned kmax = k_max;
	unsigned long long* nck;
	unsigned long long tot;
	densest* ds;
	g = getEdgeList(adj, selected);

	mkgraph(g);
	//kcore(g, kmax);
	getrankloop(g, kmax);
	relabelnodes(g);
	mkspecial(g);

	std::vector<std::vector<int>>clique_out;
	onepassfordoubletri(g, kmax, clique_out, v2cliques);

	free(g->edges);
	free(g->rank);
	free(g->d);
	free(g->cd);
	free(g->adj);
	free(g->d0);
	free(g->cd0);
	free(g->adj0);
	free(g);
	return clique_out;
}

std::vector<std::vector<int>>get_tritriangle(std::vector<std::vector<int>> adj, std::vector<std::vector<int>>& v2cliques, unsigned k_max, std::vector<bool>selected)
{
	sparse* g;
	unsigned kmax = k_max;
	unsigned long long* nck;
	unsigned long long tot;
	densest* ds;
	g = getEdgeList(adj, selected);

	mkgraph(g);
	//kcore(g, kmax);
	getrankloop(g, kmax);
	relabelnodes(g);
	mkspecial(g);

	std::vector<std::vector<int>>clique_out;
	onepassfortritri(g, kmax, clique_out, v2cliques);

	free(g->edges);
	free(g->rank);
	free(g->d);
	free(g->cd);
	free(g->adj);
	free(g->d0);
	free(g->cd0);
	free(g->adj0);
	free(g);
	return clique_out;
}

