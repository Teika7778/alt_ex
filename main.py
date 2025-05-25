import networkx as nx
import matplotlib.pyplot as plt
import copy

h = {}
index = {}
G = {}

# отображение из множества строк нуклеотидов в множество
# ребер сжатого графа. Вершина -> набор ребер

k = 14


class Edge:
    def __init__(self, seq):
        self.seq = seq
        self.cov = 0

    def __str__(self):
        return str(self.seq)


def get_all_out(a: str) -> list:
    return [key for key in index if len(key) == k and key[:-1] == a]


def get_all_in(a: str) -> list:
    return [key for key in index if len(key) == k and key[1:] == a]


def get_substrings(s, k):
    if k <= 0 or k > len(s):
        return []
    return [s[i:i + k] for i in range(len(s) - k + 1)]


def proc_str(str):
    p = get_substrings(str, k)
    for s in p:
        rev = reverse_complement(s)
        if s not in index:
            index[s] = None
            index[rev] = None


def reverse_complement(s):
    reversed_str = s[::-1]

    complement = {
        'T': 'A',
        'G': 'C',
        'A': 'T',
        'C': 'G'
    }
    return ''.join([complement.get(char, char) for char in reversed_str])


def readfile():
    with open('ebola_read.fq', 'r') as file:
        for i, line in enumerate(file):
            if i % 4 != 1:
                continue
            line = line.replace("\n", "")
            proc_str(line)


def count_cov():
    with open('ebola_read.fq', 'r') as file:
        for i, line in enumerate(file):
            if i % 4 != 1:
                continue
            line = line.replace("\n", "")
            proc_cov(line)


def proc_cov(st):
    p = get_substrings(st, k)
    for s in p:
        rev = reverse_complement(s)
        index[s][0].cov += 1 / len(index[s][0].seq)
        index[rev][0].cov += 1 / len(index[rev][0].seq)


def build_graph():
    for u in index:
        if index[u] is not None:
            continue
        # имеем текущее РЕБРО u, попробуем перейти на то ребро которое справа
        w = u
        right_v = u[:-1]  # найдем единственную вершину которая может стоять справа
        # следующее ребро нужно добавлять ТОЛЬКО когда у вершины него есть ровно 1 входящее и 1 исходящее
        while len(get_all_in(right_v)) == 1 and len(get_all_out(right_v)) == 1:
            w = get_all_in(right_v)[0]
            right_v = get_all_in(right_v)[0][:-1]
        # таким образом в конце цикла у нас будет вершина которая должна быть начальной для ребра сжатого графа
        left_v = w[1:]  # то, куда ведет наше ребро - это
        path = [w]
        while len(get_all_in(left_v)) == 1 and len(get_all_out(left_v)) == 1:
            w = get_all_out(left_v)[0]
            path.append(w)
            left_v = get_all_out(left_v)[0][1:]

        j = Edge(path)
        if right_v in G:
            G[right_v].append(j)
        else:
            G[right_v] = [j]
        for i, p in enumerate(path):
            index[p] = (j, i)


def edge_dest(e: Edge):
    return e.seq[-1][1:]


def edge_source(e: Edge):
    return e.seq[0][:-1]


# попытка быстро найти все входящие ребра
def all_in(v: str):
    verts = get_all_in(v)
    ret = []
    for u in verts:
        edge = index[u][0]
        s = edge_source(edge)
        if edge in G[s]:
            ret += [edge]
    return ret


def comp_edge(e):
    # получение комплиментарного ребра, для этого нужно найти вершину комплиментарную к
    # вершине в которую идет дуга
    v = reverse_complement(edge_dest(e))
    # затем достаточно перебрать только первые элементы последовательностей ребер
    # так как каждая из них в графе может присутствовать только в одном месте
    key = reverse_complement(e.seq[0])
    if v not in G:
        return create_comp_edge(e)
    for u in G[v]:
        if u.seq[-1] == key:
            return u
    return create_comp_edge(e)


def create_comp_edge(e):
    # иногда оказывается что комплиметарное ребро отсутвует в графе, тогда его приедстя создавать
    seq = []
    for x in e.seq:
        seq = [reverse_complement(x)] + seq
    E = Edge(seq)
    E.cov = e.cov
    return E


def rm_edge(E):
    # удаляем ребро, а также его сопряженное ребро, для этого нужно найти вершину,
    # из которой наше ребро идет, а затем удалить объект ребра из списка всех ребер
    to_rem = {E, comp_edge(E)}
    for u in to_rem:
        if u not in G[edge_source(u)]:
            continue
        G[edge_source(u)].remove(u)


def rm_vert(V):
    # удаление вершины из графа, при этом удаляются все инц. ребра (я надеюсь)
    verts = {V, reverse_complement(V)}
    for v in verts:
        # иногда получается так, что саму вершину уже удалили
        if v in G:
            del G[v]
        # перебор всего входящих в вершину ребер
        for u in G:
            if G[u] is None:
                continue
            for e in G[u]:
                if edge_dest(e) == v:
                    G[u].remove(e)


def add_edge(E):
    # умеренно загадочная функция, исходя из последовательности в ребре легко установить
    # начальную вершину, просто запишем в нее новое ребро, при этом интересно, что этой веришны
    eds = {E, comp_edge(E)}
    for e in eds:
        if edge_source(e) in G:
            G[edge_source(e)].append(e)
        else:
            G[edge_source(e)] = [e]


def simpl(Edges):
    # заменяет последовательность ребер одним ребром
    s = []
    c = 0
    l = 0
    for e in Edges:
        c += e.cov * len(e.seq)
        l += len(e.seq)
        for n in e.seq:
            s.append(n)

    E = Edge(s)
    E.cov = c / l
    add_edge(E)
    for e in Edges:
        rm_edge(e)


def compress(v):
    if v not in G:
        return
    o = G[v]
    i = all_in(v)
    if len(o) == 1 and len(i) == 1:
        path = [i[0], o[0]]
        simpl(path)


def get_all_edges():
    e = []
    for v in G:
        for u in G[v]:
            e.append(u)
    return e


def split_array_by_ratio(arr, ratio):
    sum_ratio = sum(ratio)
    n = len(arr)
    fractions = [(r / sum_ratio) * n for r in ratio]
    sizes = [int(f) for f in fractions]
    remainder = n - sum(sizes)

    if remainder > 0:
        fractional_parts = [(i, fractions[i] - sizes[i]) for i in range(len(ratio))]
        fractional_parts.sort(key=lambda x: (-x[1], x[0]))
        for i in range(remainder):
            idx = fractional_parts[i][0]
            sizes[idx] += 1

    result = []
    current = 0
    for size in sizes:
        result.append(arr[current:current + size])
        current += size
    return result


def separate(e: Edge, a: list):
    c = e.cov
    arr = e.seq
    s = split_array_by_ratio(arr, a)
    rm_edge(e)
    for t in s:
        p = Edge(t)
        p.cov = c
        add_edge(p)


def rm_buble(r, e_cov, lim=3):
    P = get_all_edges()
    while len(P) != 0:
        p = P[0]
        start = edge_source(p)
        end = edge_dest(p)
        # запускаем небольшой поиск в ширину из вершины начала
        pathes = [[start]]
        while all(x[-1] != end for x in pathes) and all(len(x) < lim for x in pathes):
            np = []
            for pa in pathes:
                for e in G[pa[-1]]:
                    d = edge_dest(e)
                    updated = copy.deepcopy(pa)
                    updated.append(d)
                    np.append(updated)
            pathes += np
            for x in pathes:
                if x[-1] == end:
                    # получается, что найден пузырь
                    pass
        P.remove(p)


def rm_tips(l_max, r, e_cov):
    P = get_all_edges()
    while len(P) != 0:
        p = P[0]
        if p not in G[edge_source(p)]:
            P.remove(p)
            continue
        if len(p.seq) >= l_max:
            P.remove(p)
            continue
        d = edge_dest(p)
        s = edge_source(p)
        if d in G and len(G[d]) != 0:
            P.remove(p)
            continue
        w = None
        for e in G[s]:
            if e == p:
                continue
            if (w is None) or (w.cov < e.cov):
                w = e
        if w is None:
            P.remove(p)
            continue
        if not (p.cov > e_cov or p.cov * r > w.cov):
            rm_edge(p)

        P.remove(p)


def rm_chem(cov):
    P = get_all_edges()
    while len(P) != 0:
        p = P[0]
        if p.cov < cov:
            rm_edge(p)
            compress(edge_dest(p))
            compress(edge_source(p))

        P.remove(p)


def draw_compressed_graph(adjacency):
    G = nx.MultiDiGraph()
    for source, edges in adjacency.items():
        if edges is None:
            continue
        for edge in edges:
            if not edge.seq:
                continue
            target = edge.seq[-1][1:]  # Последний символ -> целевая вершина
            G.add_edge(source, target, label=str(edge.cov))

    pos = nx.spring_layout(G, seed=22)
    nx.draw(
        G, pos,
        with_labels=True,
        node_size=2500,
        node_color="lightgreen",
        font_size=12,
        arrows=True,
        connectionstyle="arc3,rad=0.2"
    )

    edge_labels = {(u, v, k): d['label']
                   for u, v, k, d in G.edges(keys=True, data=True)}
    nx.draw_networkx_edge_labels(
        G, pos,
        edge_labels=edge_labels,
        font_color="darkred",
        label_pos=0.75
    )

    plt.title("Сжатый граф де Брюина")
    plt.show()


readfile()
build_graph()
draw_compressed_graph(G)
count_cov()
rm_chem(2)
rm_tips(150, 2, 30)
rm_buble(2, 30)
draw_compressed_graph(G)

for source, edges in G.items():
    if edges is None:
        continue
    for edge in edges:
        if not edge.seq:
            continue
        if len(get_all_in(source)) == 0:
            s = source
            cur = source
            while cur in G and len(G[cur]) == 1:
                for n in G[cur][0].seq:
                    s += n[-1]
                cur = edge_dest(G[cur][0])
            print(s)

