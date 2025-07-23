
from triangular import *


# tested 1..30
def test_young_2up(n):
    for tp in TriangularPartitions(n):
        if not len(list(tp.up())) <= 2:
            return tp
    return True

def test_symmetry_distribution(d):
    return all(d.get((y,x),0) == d.get((x,y),0) for x,y in d)



def test_all_partitions_symmetry_distance_sim_intervals(n):
    for p in TriangularPartitions(n):
        print(p)
        if not test_symmetry_distribution(p.distance_sim_distribution_intervals()):
            return p
    return True

def test_deficit_cells(tdp):
    return all(tdp.is_sim_cell(i,j) != tdp.is_deficit_cell(i,j) for i,j in tdp.path().cells())

def test_deficit_cells_partition(p):
    for tdp in p.triangular_dyck_paths():
        if not test_deficit_cells(tdp):
            return tdp
    return True

# tested 3..12
def test_all_partitions_deficit_cells(n):
    for p in TriangularPartitions(n):
        print(p)
        r = test_deficit_cells_partition(p)
        if r != True:
            return r
    return True

# tested 3..21
# def test_all_partitions_natural_down(n):
    # for p in TriangularPartitions(n):
        # print(p)
        # if p.natural_down() is None:
            # return p
    # return True

# No for 5,3,2
def test_all_partitions_diagonal_orientation(n):
    lines = {TriangularPartition([n]), TriangularPartition([1]*n)}
    for p in TriangularPartitions(n):
        print(p)
        if not p in lines and p.diagonal_orientation() is None:
            return p
    return True

# False for 12
# 5,4,2,1
# 4,3,2,2,1
def test_all_partitions_diagonal_oriented(n):
    return all(t.is_diagonal_oriented() for t in TriangularPartitions(n))


def test_vertical_horizontal_deficit(tdp):
    return all(tdp.is_horizontal_deficit(i,j) != tdp.is_vertical_deficit(i,j) for i,j in tdp.deficit_cells())

def test_vertical_horizontal_deficit_partition(p):
    for tdp in p.triangular_dyck_paths():
        if not test_vertical_horizontal_deficit(tdp):
            return tdp
    return True

# tested 3..12
def test_all_partitions_vertical_horizontal_deficit(n):
    for p in TriangularPartitions(n):
        print(p)
        r = test_vertical_horizontal_deficit_partition(p)
        if r != True:
            return r
    return True

def test_corner_vertical_horizontal_deficit(tdp):
    for i,j in tdp.corners():
        if i == 0 or j == 0:
            return True
        if tdp.is_horizontal_deficit(i,j-1) and tdp.is_vertical_deficit(i-1,j):
            return False
    return True

def test_corner_vertical_horizontal_deficit_partition(p):
    for tdp in p.triangular_dyck_paths():
        if not test_corner_vertical_horizontal_deficit(tdp):
            return tdp
    return True

# tested 3..12
def test_all_partitions_corner_vertical_horizontal_deficit(n):
    for p in TriangularPartitions(n):
        print(p)
        r = test_corner_vertical_horizontal_deficit_partition(p)
        if r != True:
            return r
    return True

# broken 5,3,1 fixed
# broken 6,3,1 fixed
# broken 6,4,2,1 fixed
# tested 3..18
# OLD VERSION -- NOT SYMMETRIC
def test_all_lattice_trianglar_tamari(n):
    for p in TriangularPartitions(n):
        print(p)
        if not p.triangular_tamari_poset().is_lattice():
            return p
    return True

def test_all_lattice_immediate_deficit_tamari(n):
    for p in TriangularPartitions(n):
        print(p)
        if not p.immediate_deficit_tamari_poset().is_lattice():
            return p
    return True


def ideal_poset(P,x,y):
    return Poset(([v for v in P if P.le(v,x) and P.le(v,y)],lambda a,b: P.le(a,b)))

# tested 1..9
# Error for 6,3,1
def test_all_lattice_unambiguous_trianglar_tamari(n):
    for p in TriangularPartitions(n):
        if p.is_unambiguous():
            print(p)
            if not p.triangular_tamari_poset().is_lattice():
                return p
    return True

def test_all_distribution_unambiguous_triangular_tamari(n):
    for p in TriangularPartitions(n):
        if p.is_unambiguous():
            print(p)
            P = p.triangular_tamari_poset()
            d = poset_distance_sim_distribution_intervals(P)
            if not test_symmetry_distribution(d):
                return p
    return True

def face_2d(L):
    for e in L:
        U = L.upper_covers(e)
        for i in range(len(U)-1):
            for j in range(i+1,len(U)):
                yield L.interval(e, L.join(U[i],U[j]))

def is_pentagonal(L):
    return all(len(x) == 5 or len(x) == 4 for x in face_2d(L))

def poset_distance(P, dw1, dw2):
    I = P.subposet(P.interval(dw1, dw2))
    return max(len(c) for c in I.maximal_chains()) - 1

def poset_distance_sim_distribution_intervals(P):
    d = {}
    for dp1, dp2 in P.relations():
        distance = poset_distance(P, dp1, dp2)
        sim = dp2.sim()
        d[(distance, sim)] = d.get((distance, sim),0) + 1
    return d

def intervals_polynomial(P, q, t):
    d = poset_distance_sim_distribution_intervals(P)
    return sum(d[k] * q**k[0] * t**k[1] for k in d)

def schur_poset_distance_sim_area_distribution_interval(P):
    d = {}
    for dp1, dp2 in P.relations():
        I = P.subposet(P.interval(dp1,dp2))
        #if dp1.deficit() > 0 and dp2.deficit() > 0 and any(dp.deficit() == 0 for dp in I):
            #continue
        distance = max(len(c) for c in I.maximal_chains()) - 1
        sim = dp2.sim()
        area = dp1.area()
        d[(distance, sim,area)] = d.get((distance, sim,area),0) + 1
    return d

def schur_poset_polynomial(P, q, t, r):
    d = schur_poset_distance_sim_area_distribution_interval(P)
    return sum(d[k] * q**k[0] * t**k[1] * r**k[2] for k in d)

# tested up to size 20
def test_corners_maximal(tr):
    L = tr.triangular_standard_labels()
    corners = set(tr.corners())
    not_corners = set(c for c in tr.cells() if not c in corners)
    return all(L[c1] < L[c2] for c1 in not_corners for c2 in corners)


### r param

# def r_param(P, dp1, dp2):
    # C1 = set(dp1.lead_deficit_cells())
    # C2 = set(dp2.lead_deficit_cells())
    # return dp2.area() - poset_distance(P,dp1,dp2) - sum(1 for c1 in C1 if not c1 in C2)

# def equilibrum(P,dp):
    # return dp.sim() - poset_distance(P,dp,P.maximal_elements()[0])

# def r_param(P,dp1,dp2):
    # eq1 = equilibrum(P,dp1)
    # eq2 = equilibrum(P,dp2)
    # return dp1.area() + ((eq1 - eq2) if abs(eq2) < abs(eq1) else 0)

# def r_param(dp1,dp2):
    # rcells1 = set(dp1.triangular_tamari_removable_cells())
    # plus_cells1 = set(c for c in dp1.similar_cells() if not c in rcells1)
    # minus_cells1 = set(c for c in rcells1 if dp1.is_deficit_cell(*c))


    # rcells2 = set(dp2.triangular_tamari_removable_cells())
    # plus_cells2 = set(c for c in dp2.similar_cells() if not c in rcells2)
    # minus_cells2 = set(c for c in rcells2 if dp2.is_deficit_cell(*c))

    # return dp1.area() + len(plus_cells1.difference(plus_cells2)) - len(minus_cells1.difference(minus_cells2))

def r_param(dp1, dp2):
    return dp1.area() + len(dp1.triangular_tamari_plus_cells().difference(dp2.triangular_tamari_plus_cells())) - len(dp1.triangular_tamari_minus_cells().difference(dp2.triangular_tamari_minus_cells()))

def q_t_r_tamari_distribution(P):
    d = {}
    for dp1, dp2 in P.relations():
        sim = dp2.sim()
        distance = poset_distance(P,dp1,dp2)
        rr = r_param(dp1,dp2)
        d[(distance,sim,rr)] = d.get((distance,sim,rr),0) + 1
    return d

def q_t_r_tamari_polynomal(P,q,t,r):
    d = q_t_r_tamari_distribution(P)
    return sum(d[k] * q**k[0]*t**k[1]*r**k[2] for k in d)

## rational case

def test_symmetry_rational_partitions(n):
    for tr in TriangularPartitions(n):
        if tr.length() > 1 and tr[0] >= tr.length() and tr.is_rational():
            print(tr)
            d = tr.path_distance_sim_distribution_intervals()
            if not test_symmetry_distribution(d):
                return tr
    return True

## 2 lines

def simple_q_t_r_distribution(P):
    d = {}
    for dp1, dp2 in P.relations():
        sim = dp2.sim()
        distance = poset_distance(P,dp1,dp2)
        rr = dp1.area()
        d[(distance,sim,rr)] = d.get((distance,sim,rr),0) + 1
    return d

def simple_q_t_r_polynomal(P,q,t,r):
    d = simple_q_t_r_distribution(P)
    return sum(d[k] * q**k[0]*t**k[1]*r**k[2] for k in d)

# tested 4 -- 13
def test_2_lines_triangular_tamari(n):
    for tr in TriangularPartitions(n, length=2):
        d = tr.triangular_distance_sim_distribution_intervals()
        if not test_symmetry_distribution(d):
            print(tr)
            return False
    return True

def test_2_lines_path_tamari(n):
    for tr in TriangularPartitions(n, length=2):
        d = tr.path_distance_sim_distribution_intervals()
        if not test_symmetry_distribution(d):
            print(tr)
            return False
    return True


## top down Labels

def top_down_q_t_distribution(tr, tableau = None):
    if tableau is None:
        tableau = tr.top_down_standard_tableau()
    d = {}
    for dp in tr.triangular_dyck_paths():
        area = dp.area()
        sim = tr.size() - area - dp.deficit(tableau = tableau)
        d[(area,sim)] = d.get((area,sim),0) + 1
    return d

def PRV_poset_distribution(tr, P, tableau = None):
    if tableau is None:
        tableau = tr.top_down_standard_tableau()
    labels = {c: tableau.entry(c) for c in tr.cells()}
    d = {}
    for dp1, dp2 in P.relations():
        sim = tr.size() - dp2.area() - dp2.deficit(tableau = tableau)
        distance = poset_distance(P,dp1,dp2)
        d[(distance,sim)] = d.get((distance,sim),0) + 1
    return d

def symmetric_labelings(tr):
    for t in tr.partition().standard_tableaux():
        d = tr.area_sim_distribution(t)
        if test_symmetry_distribution(d):
            yield t

# tested 3 .. 21
def test_both_symmetry(n):
    for tr in TriangularPartitions(n):
        d = tr.top_down_area_sim_distribution()
        if test_symmetry_distribution(d):
            print(tr)
            P = tr.path_tamari_lattice()
            d2 = PRV_poset_distribution(tr,P)
            if not test_symmetry_distribution(d2):
                return tr
    return True

# tested 3..21
def test_schur_positive_PRVsimsym(n, q, t):
    for tr in TriangularPartitions(n):
        d = PRV_q_t_distribution(tr)
        if test_symmetry_distribution(d):
            print(tr)
            P = tr.path_tamari_lattice()
            d2 = PRV_poset_distribution(tr,P)
            pol = sum(d2[k]*q**k[0]*t**k[1] for k in d2)
            S = schur.from_polynomial(pol)
            if not all(c>0 for p,c in S if len(P) <= 2):
                return tr
    return True


# tested 3..22
def test_slope_labels(n):
    for tr in TriangularPartitions(n):
        for labels in tr.slope_labels_distincts():
            d = PRV_q_t_distribution(tr, labels)
            if not test_symmetry_distribution(d):
                return tr
    return True

### sim sym

def q_t_distribution(tr, tableau = None):
    labels = {c: tableau.entry(c) for c in tr.cells()}
    d = {}
    for dp in tr.triangular_dyck_paths():
        area = dp.area()
        sim = tr.size() - area - dp.deficit(labels = labels)
        d[(area,sim)] = d.get((area,sim),0) + 1
    return d

def k2simsym(tr):
    m,n = tr[0],tr[1]

    for k in range(m - 2*n +2):
        row2 = [m+n-k - 2*i for i in range(n)]
        row2.reverse()
        row1 = [i for i in range(1,m+n+1) if not i in set(row2)]
        yield Tableau([row1,row2])

    if n == 2:
        row2 = [2,5]
        row1 = [i for i in range(1,m+n+1) if not i in set(row2)]
        yield Tableau([row1,row2])

# tested 3..13
def test_k2simsym(n):
    for tr in TriangularPartitions(n, length = 2):
        S1 = set(k2simsym(tr))
        S2 = set(symmetric_labelings(tr))
        if S1 != S2:
            return tr
    return True


def test_PRV_polynomials(pols, q, t):
    wrong = []
    for key in pols:
        tr = TriangularPartition(key)
        d = PRV_q_t_distribution(tr)
        if test_symmetry_distribution(d):
            print(tr)
            P = tr.path_tamari_lattice()
            d2 = PRV_poset_distribution(tr,P)
            pol1 = sum(d2[k]*q**k[0]*t**k[1] for k in d2)
            pol2 = pols[key].expand(3)(q,t,1)
            if not pol1 == pol2:
                wrong.append(tr)
    return wrong

Sym = SymmetricFunctions(QQ)
schur = Sym.schur()

SCHURS = {
(0,) : schur[0],
(1,) : schur[1],
(2,) : schur[2],
(2, 1) : schur[1, 1] + schur[3],
(3,) : schur[3],
(3, 1) : schur[2, 1] + schur[4],
(4,) : schur[4],
(3, 2) : schur[3, 1] + schur[5],
(4, 1) : schur[3, 1] + schur[5],
(5,) : schur[5],
(3, 2, 1) : schur[1, 1, 1] + schur[3, 1] + schur[4, 1] + schur[6],
(4, 2) : schur[2, 2] + schur[4, 1] + schur[6],
(5, 1) : schur[4, 1] + schur[6],
(6,) : schur[6],
(4, 2, 1) : schur[2, 1, 1] + schur[3, 2] + schur[4, 1] + schur[5, 1] + schur[7],
(5, 2) : schur[3, 2] + schur[5, 1] + schur[7],
(6, 1) : schur[5, 1] + schur[7],
(7,) : schur[7],
(4, 3, 1) : schur[3, 1, 1] + schur[4, 2] + schur[5, 1] + schur[6, 1] + schur[8],
(5, 3) : schur[4, 2] + schur[6, 1] + schur[8],
(6, 2) : schur[4, 2] + schur[6, 1] + schur[8],
(7, 1) : schur[6, 1] + schur[8],
(8,) : schur[8],
(4, 3, 2) : schur[3, 3] + schur[4, 1, 1] + schur[5, 2] + schur[6, 1] + schur[7, 1] + schur[9],
(5, 3, 1) : schur[2, 2, 1] + schur[4, 1, 1] + schur[4, 2] + schur[5, 2] + schur[6, 1] + schur[7, 1] + schur[9],
(6, 3) : schur[3, 3] + schur[5, 2] + schur[7, 1] + schur[9],
(7, 2) : schur[5, 2] + schur[7, 1] + schur[9],
(8, 1) : schur[7, 1] + schur[9],
(9,) : schur[9],
(4, 3, 2, 1) : schur[1, 1, 1, 1] + schur[3, 1, 1] + schur[4, 1, 1] + schur[4, 2] + schur[4, 3] + schur[5, 1, 1] + schur[6, 1] + schur[6, 2] + schur[7, 1] + schur[8, 1] + schur[10],
(5, 3, 2) : schur[3, 2, 1] + schur[4, 3] + schur[5, 1, 1] + schur[5, 2] + schur[6, 2] + schur[7, 1] + schur[8, 1] + schur[10],
(6, 3, 1) : schur[3, 2, 1] + schur[4, 3] + schur[5, 1, 1] + schur[5, 2] + schur[6, 2] + schur[7, 1] + schur[8, 1] + schur[10],
(7, 3) : schur[4, 3] + schur[6, 2] + schur[8, 1] + schur[10],
(8, 2) : schur[6, 2] + schur[8, 1] + schur[10],
(9, 1) : schur[8, 1] + schur[10],
(10,) : schur[10],
(5, 3, 2, 1) : schur[2, 1, 1, 1] + schur[3, 2, 1] + schur[4, 1, 1] + schur[4, 2, 1] + schur[4, 3] + schur[5, 1, 1] + schur[5, 2] + schur[5, 3] + schur[6, 1, 1] + schur[6, 2] + schur[7, 1] + schur[7, 2] + schur[8, 1] + schur[9, 1] + schur[11],
(5, 4, 2) : schur[4, 2, 1] + schur[5, 3] + schur[6, 1, 1] + schur[6, 2] + schur[7, 2] + schur[8, 1] + schur[9, 1] + schur[11],
(6, 4, 1) : schur[4, 2, 1] + schur[5, 3] + schur[6, 1, 1] + schur[6, 2] + schur[7, 2] + schur[8, 1] + schur[9, 1] + schur[11],
(7, 4) : schur[5, 3] + schur[7, 2] + schur[9, 1] + schur[11],
(8, 3) : schur[5, 3] + schur[7, 2] + schur[9, 1] + schur[11],
(9, 2) : schur[7, 2] + schur[9, 1] + schur[11],
(10, 1) : schur[9, 1] + schur[11],
(11,) : schur[11],
(5, 4, 2, 1) : schur[3, 1, 1, 1] + schur[3, 3, 1] + schur[4, 2, 1] + schur[4, 4] + schur[5, 1, 1] + schur[5, 2, 1] + schur[5, 3] + schur[6, 1, 1] + schur[6, 2] + schur[6, 3] + schur[7, 1, 1] + schur[7, 2] + schur[8, 1] + schur[8, 2] + schur[9, 1] + schur[10, 1] + schur[12],
(6, 4, 2) : schur[2, 2, 2] + schur[4, 2, 1] + schur[4, 4] + schur[5, 2, 1] + schur[6, 2] + schur[6, 3] + schur[7, 1, 1] + schur[7, 2] + schur[8, 2] + schur[9, 1] + schur[10, 1] + schur[12],
(7, 4, 1) : schur[3, 3, 1] + schur[5, 2, 1] + schur[5, 3] + schur[6, 3] + schur[7, 1, 1] + schur[7, 2] + schur[8, 2] + schur[9, 1] + schur[10, 1] + schur[12],
(8, 4) : schur[4, 4] + schur[6, 3] + schur[8, 2] + schur[10, 1] + schur[12],
(9, 3) : schur[6, 3] + schur[8, 2] + schur[10, 1] + schur[12],
(10, 2) : schur[8, 2] + schur[10, 1] + schur[12],
(11, 1) : schur[10, 1] + schur[12],
(12,) : schur[12],
(5, 4, 3, 1) : schur[4, 1, 1, 1] + schur[4, 3, 1] + schur[5, 2, 1] + schur[5, 4] + schur[6, 1, 1] + schur[6, 2, 1] + schur[6, 3] + schur[7, 1, 1] + schur[7, 2] + schur[7, 3] + schur[8, 1, 1] + schur[8, 2] + schur[9, 1] + schur[9, 2] + schur[10, 1] + schur[11, 1] + schur[13],
(6, 4, 2, 1) : schur[2, 2, 1, 1] + schur[4, 1, 1, 1] + schur[4, 2, 1] + schur[4, 3, 1] + 2*schur[5, 2, 1] + schur[5, 3] + schur[5, 4] + schur[6, 1, 1] + schur[6, 2, 1] + schur[6, 3] + schur[7, 1, 1] + 2*schur[7, 2] + schur[7, 3] + schur[8, 1, 1] + schur[8, 2] + schur[9, 1] + schur[9, 2] + schur[10, 1] + schur[11, 1] + schur[13] + schur[3,2,2],
(7, 4, 2) : schur[3, 2, 2] + schur[4, 3, 1] + schur[5, 2, 1] + schur[5, 4] + schur[6, 2, 1] + schur[6, 3] + schur[7, 2] + schur[7, 3] + schur[8, 1, 1] + schur[8, 2] + schur[9, 2] + schur[10, 1] + schur[11, 1] + schur[13],
(8, 4, 1) : schur[4, 3, 1] + schur[5, 4] + schur[6, 2, 1] + schur[6, 3] + schur[7, 3] + schur[8, 1, 1] + schur[8, 2] + schur[9, 2] + schur[10, 1] + schur[11, 1] + schur[13],
(9, 4) : schur[5, 4] + schur[7, 3] + schur[9, 2] + schur[11, 1] + schur[13],
(10, 3) : schur[7, 3] + schur[9, 2] + schur[11, 1] + schur[13],
(11, 2) : schur[9, 2] + schur[11, 1] + schur[13],
(12, 1) : schur[11, 1] + schur[13],
(13,) : schur[13],
(5, 4, 3, 2) : schur[4, 3, 1] + schur[5, 1, 1, 1] + schur[5, 3, 1] + schur[6, 2, 1] + schur[6, 3] + schur[6, 4] + schur[7, 1, 1] + schur[7, 2, 1] + schur[7, 3] + schur[8, 1, 1] + schur[8, 2] + schur[8, 3] + schur[9, 1, 1] + schur[9, 2] + schur[10, 1] + schur[10, 2] + schur[11, 1] + schur[12, 1] + schur[14],
(6, 4, 3, 1) : schur[3, 2, 1, 1] + schur[4, 3, 1] + schur[5, 1, 1, 1] + schur[5, 2, 1] + schur[5, 3, 1] + schur[5, 4] + 2*schur[6, 2, 1] + schur[6, 3] + schur[6, 4] + schur[7, 1, 1] + schur[7, 2, 1] + schur[7, 3] + schur[8, 1, 1] + 2*schur[8, 2] + schur[8, 3] + schur[9, 1, 1] + schur[9, 2] + schur[10, 1] + schur[10, 2] + schur[11, 1] + schur[12, 1] + schur[14],
(7, 5, 2) : schur[4, 2, 2] + schur[5, 3, 1] + schur[6, 2, 1] + schur[6, 4] + schur[7, 2, 1] + schur[7, 3] + schur[8, 2] + schur[8, 3] + schur[9, 1, 1] + schur[9, 2] + schur[10, 2] + schur[11, 1] + schur[12, 1] + schur[14],
(8, 5, 1) : schur[5, 3, 1] + schur[6, 4] + schur[7, 2, 1] + schur[7, 3] + schur[8, 3] + schur[9, 1, 1] + schur[9, 2] + schur[10, 2] + schur[11, 1] + schur[12, 1] + schur[14],
(9, 5) : schur[6, 4] + schur[8, 3] + schur[10, 2] + schur[12, 1] + schur[14],
(10, 4) : schur[6, 4] + schur[8, 3] + schur[10, 2] + schur[12, 1] + schur[14],
(11, 3) : schur[8, 3] + schur[10, 2] + schur[12, 1] + schur[14],
(12, 2) : schur[10, 2] + schur[12, 1] + schur[14],
(13, 1) : schur[12, 1] + schur[14],
(14,) : schur[14],
(5, 4, 3, 2, 1) : schur[1, 1, 1, 1, 1] + schur[3, 1, 1, 1] + schur[4, 1, 1, 1] + schur[4, 2, 1] + schur[4, 3, 1] + schur[4, 4] + schur[4, 4, 1] + schur[5, 1, 1, 1] + schur[5, 2, 1] + schur[5, 3, 1] + schur[6, 1, 1] + schur[6, 1, 1, 1] + schur[6, 2, 1] + schur[6, 3] + schur[6, 3, 1] + schur[6, 4] + schur[7, 1, 1] + schur[7, 2] + schur[7, 2, 1] + schur[7, 3] + schur[7, 4] + 2*schur[8, 1, 1] + schur[8, 2] + schur[8, 2, 1] + schur[8, 3] + schur[9, 1, 1] + schur[9, 2] + schur[9, 3] + schur[10, 1] + schur[10, 1, 1] + schur[10, 2] + schur[11, 1] + schur[11, 2] + schur[12, 1] + schur[13, 1] + schur[15],
(6, 4, 3, 2) : schur[4, 2, 1, 1] + schur[4, 4, 1] + schur[5, 3, 1] + schur[5, 4, 1] + schur[5, 5] + schur[6, 1, 1, 1] + schur[6, 2, 1] + 2*schur[6, 3, 1] + schur[6, 4] + 2*schur[7, 2, 1] + 2*schur[7, 3] + schur[7, 4] + schur[8, 1, 1] + schur[8, 2, 1] + schur[8, 3] + schur[9, 1, 1] + 2*schur[9, 2] + schur[9, 3] + schur[10, 1, 1] + schur[10, 2] + schur[11, 1] + schur[11, 2] + schur[12, 1] + schur[13, 1] + schur[15],
(6, 5, 3, 1) : schur[4, 2, 1, 1] + schur[4, 4, 1] + schur[5, 3, 1] + schur[5, 5] + schur[6, 1, 1, 1] + schur[6, 2, 1] + 2*schur[6, 3, 1] + schur[6, 4] + 2*schur[7, 2, 1] + schur[7, 3] + schur[7, 4] + schur[8, 1, 1] + schur[8, 2, 1] + schur[8, 3] + schur[9, 1, 1] + 2*schur[9, 2] + schur[9, 3] + schur[10, 1, 1] + schur[10, 2] + schur[11, 1] + schur[11, 2] + schur[12, 1] + schur[13, 1] + schur[15],
(7, 5, 3) : schur[4, 4, 1] + schur[5, 2, 2] + schur[6, 3, 1] + schur[6, 4] + schur[7, 2, 1] + schur[7, 4] + schur[8, 2, 1] + schur[8, 3] + schur[9, 2] + schur[9, 3] + schur[10, 1, 1] + schur[10, 2] + schur[11, 2] + schur[12, 1] + schur[13, 1] + schur[15],
(8, 5, 2) : schur[3, 3, 2] + schur[5, 2, 2] + schur[5, 3, 1] + schur[5, 5] + schur[6, 3, 1] + schur[7, 2, 1] + schur[7, 3] + schur[7, 4] + schur[8, 2, 1] + schur[8, 3] + schur[9, 2] + schur[9, 3] + schur[10, 1, 1] + schur[10, 2] + schur[11, 2] + schur[12, 1] + schur[13, 1] + schur[15],
(9, 5, 1) : schur[4, 4, 1] + schur[6, 3, 1] + schur[6, 4] + schur[7, 4] + schur[8, 2, 1] + schur[8, 3] + schur[9, 3] + schur[10, 1, 1] + schur[10, 2] + schur[11, 2] + schur[12, 1] + schur[13, 1] + schur[15],
(10, 5) : schur[5, 5] + schur[7, 4] + schur[9, 3] + schur[11, 2] + schur[13, 1] + schur[15],
(11, 4) : schur[7, 4] + schur[9, 3] + schur[11, 2] + schur[13, 1] + schur[15],
(12, 3) : schur[9, 3] + schur[11, 2] + schur[13, 1] + schur[15],
(13, 2) : schur[11, 2] + schur[13, 1] + schur[15],
(14, 1) : schur[13, 1] + schur[15],
(15,) : schur[15],
(6, 4, 3, 2, 1) : schur[2, 1, 1, 1, 1] + schur[3, 2, 1, 1] + schur[4, 1, 1, 1] + schur[4, 2, 1, 1] + schur[4, 2, 2] + schur[4, 3, 1] + schur[4, 3, 2] + schur[4, 4, 1] + schur[5, 1, 1, 1] + schur[5, 2, 1] + schur[5, 2, 1, 1] + 2*schur[5, 3, 1] + schur[5, 4] + schur[5, 4, 1] + schur[6, 1, 1, 1] + 2*schur[6, 2, 1] + schur[6, 2, 2] + 2*schur[6, 3, 1] + schur[6, 4] + schur[6, 5] + schur[7, 1, 1] + schur[7, 1, 1, 1] + 2*schur[7, 2, 1] + 2*schur[7, 3] + schur[7, 3, 1] + schur[7, 4] + schur[8, 1, 1] + schur[8, 2] + 2*schur[8, 2, 1] + 2*schur[8, 3] + schur[8, 4] + 2*schur[9, 1, 1] + schur[9, 2] + schur[9, 2, 1] + schur[9, 3] + schur[10, 1, 1] + 2*schur[10, 2] + schur[10, 3] + schur[11, 1] + schur[11, 1, 1] + schur[11, 2] + schur[12, 1] + schur[12, 2] + schur[13, 1] + schur[14, 1] + schur[16],
(6, 5, 3, 2) : schur[3, 3, 1, 1] + schur[5, 2, 1, 1] + schur[5, 3, 1] + schur[5, 4, 1] + 2*schur[6, 3, 1] + schur[6, 4] + schur[6, 5] + schur[7, 1, 1, 1] + schur[7, 2, 1] + 2*schur[7, 3, 1] + schur[7, 4] + 2*schur[8, 2, 1] + 2*schur[8, 3] + schur[8, 4] + schur[9, 1, 1] + schur[9, 2, 1] + schur[9, 3] + schur[10, 1, 1] + 2*schur[10, 2] + schur[10, 3] + schur[11, 1, 1] + schur[11, 2] + schur[12, 1] + schur[12, 2] + schur[13, 1] + schur[14, 1] + schur[16],
(7, 5, 3, 1) : schur[2, 2, 2, 1] + schur[4, 2, 1, 1] + schur[4, 4, 1] + schur[5, 2, 1, 1] + schur[5, 3, 1] + schur[5, 4, 1] + schur[6, 2, 1] + schur[6, 3, 1] + schur[6, 4] + schur[6, 5] + schur[7, 1, 1, 1] + 2*schur[7, 2, 1] + schur[7, 3] + schur[7, 3, 1] + schur[7, 4] + 2*schur[8, 2, 1] + schur[8, 3] + schur[8, 4] + schur[9, 1, 1] + schur[9, 2] + schur[9, 2, 1] + schur[9, 3] + schur[10, 1, 1] + 2*schur[10, 2] + schur[10, 3] + schur[11, 1, 1] + schur[11, 2] + schur[12, 1] + schur[12, 2] + schur[13, 1] + schur[14, 1] + schur[16] + schur[6,2,2] + schur[5,2,2] + schur[4,2,2],
(8, 5, 3) : schur[4, 3, 2] + schur[5, 4, 1] + schur[6, 2, 2] + schur[6, 3, 1] + schur[6, 5] + schur[7, 3, 1] + schur[7, 4] + schur[8, 2, 1] + schur[8, 3] + schur[8, 4] + schur[9, 2, 1] + schur[9, 3] + schur[10, 2] + schur[10, 3] + schur[11, 1, 1] + schur[11, 2] + schur[12, 2] + schur[13, 1] + schur[14, 1] + schur[16],
(9, 5, 2) : schur[4, 3, 2] + schur[5, 4, 1] + schur[6, 2, 2] + schur[6, 3, 1] + schur[6, 5] + schur[7, 3, 1] + schur[7, 4] + schur[8, 2, 1] + schur[8, 3] + schur[8, 4] + schur[9, 2, 1] + schur[9, 3] + schur[10, 2] + schur[10, 3] + schur[11, 1, 1] + schur[11, 2] + schur[12, 2] + schur[13, 1] + schur[14, 1] + schur[16],
(10, 5, 1) : schur[5, 4, 1] + schur[6, 5] + schur[7, 3, 1] + schur[7, 4] + schur[8, 4] + schur[9, 2, 1] + schur[9, 3] + schur[10, 3] + schur[11, 1, 1] + schur[11, 2] + schur[12, 2] + schur[13, 1] + schur[14, 1] + schur[16],
(11, 5) : schur[6, 5] + schur[8, 4] + schur[10, 3] + schur[12, 2] + schur[14, 1] + schur[16],
(12, 4) : schur[8, 4] + schur[10, 3] + schur[12, 2] + schur[14, 1] + schur[16],
(13, 3) : schur[10, 3] + schur[12, 2] + schur[14, 1] + schur[16],
(14, 2) : schur[12, 2] + schur[14, 1] + schur[16],
(15, 1) : schur[14, 1] + schur[16],
(16,) : schur[16],
(6, 5, 3, 2, 1) : schur[3, 1, 1, 1, 1] + schur[3, 3, 1, 1] + schur[4, 2, 1, 1] + schur[4, 3, 1, 1] + schur[4, 3, 2] + schur[4, 4, 1] + schur[5, 1, 1, 1] + schur[5, 2, 1, 1] + schur[5, 2, 2] + schur[5, 3, 1] + schur[5, 3, 2] + 2*schur[5, 4, 1] + schur[6, 1, 1, 1] + schur[6, 2, 1] + schur[6, 2, 1, 1] + 3*schur[6, 3, 1] + schur[6, 4] + schur[6, 4, 1] + schur[6, 5] + schur[7, 1, 1, 1] + 2*schur[7, 2, 1] + schur[7, 2, 2] + 2*schur[7, 3, 1] + 2*schur[7, 4] + schur[7, 5] + schur[8, 1, 1] + schur[8, 1, 1, 1] + 2*schur[8, 2, 1] + 2*schur[8, 3] + schur[8, 3, 1] + schur[8, 4] + schur[9, 1, 1] + schur[9, 2] + 2*schur[9, 2, 1] + 2*schur[9, 3] + schur[9, 4] + 2*schur[10, 1, 1] + schur[10, 2] + schur[10, 2, 1] + schur[10, 3] + schur[11, 1, 1] + 2*schur[11, 2] + schur[11, 3] + schur[12, 1] + schur[12, 1, 1] + schur[12, 2] + schur[13, 1] + schur[13, 2] + schur[14, 1] + schur[15, 1] + schur[17],
(6, 5, 4, 2, 1) : schur[4, 1, 1, 1, 1] + schur[4, 3, 1, 1] + schur[4, 4, 2] + schur[5, 2, 1, 1] + schur[5, 3, 1, 1] + schur[5, 3, 2] + schur[5, 4, 1] + schur[5, 5, 1] + schur[6, 1, 1, 1] + schur[6, 2, 1, 1] + schur[6, 2, 2] + schur[6, 3, 1] + schur[6, 3, 2] + 2*schur[6, 4, 1] + schur[6, 6] + schur[7, 1, 1, 1] + schur[7, 2, 1] + schur[7, 2, 1, 1] + 3*schur[7, 3, 1] + schur[7, 4] + schur[7, 4, 1] + schur[7, 5] + schur[8, 1, 1, 1] + 2*schur[8, 2, 1] + schur[8, 2, 2] + 2*schur[8, 3, 1] + 2*schur[8, 4] + schur[8, 5] + schur[9, 1, 1] + schur[9, 1, 1, 1] + 2*schur[9, 2, 1] + 2*schur[9, 3] + schur[9, 3, 1] + schur[9, 4] + schur[10, 1, 1] + schur[10, 2] + 2*schur[10, 2, 1] + 2*schur[10, 3] + schur[10, 4] + 2*schur[11, 1, 1] + schur[11, 2] + schur[11, 2, 1] + schur[11, 3] + schur[12, 1, 1] + 2*schur[12, 2] + schur[12, 3] + schur[13, 1] + schur[13, 1, 1] + schur[13, 2] + schur[14, 1] + schur[14, 2] + schur[15, 1] + schur[16, 1] + schur[18],
(6, 5, 4, 3, 1) : schur[4, 3, 3] + schur[4, 4, 1, 1] + schur[5, 1, 1, 1, 1] + schur[5, 3, 1, 1] + schur[5, 4, 2] + schur[6, 2, 1, 1] + schur[6, 3, 1, 1] + schur[6, 3, 2] + 2*schur[6, 4, 1] + schur[6, 5, 1] + schur[7, 1, 1, 1] + schur[7, 2, 1, 1] + schur[7, 2, 2] + schur[7, 3, 1] + schur[7, 3, 2] + 2*schur[7, 4, 1] + schur[7, 5] + schur[7, 6] + schur[8, 1, 1, 1] + schur[8, 2, 1] + schur[8, 2, 1, 1] + 3*schur[8, 3, 1] + schur[8, 4] + schur[8, 4, 1] + schur[8, 5] + schur[9, 1, 1, 1] + 2*schur[9, 2, 1] + schur[9, 2, 2] + 2*schur[9, 3, 1] + 2*schur[9, 4] + schur[9, 5] + schur[10, 1, 1] + schur[10, 1, 1, 1] + 2*schur[10, 2, 1] + 2*schur[10, 3] + schur[10, 3, 1] + schur[10, 4] + schur[11, 1, 1] + schur[11, 2] + 2*schur[11, 2, 1] + 2*schur[11, 3] + schur[11, 4] + 2*schur[12, 1, 1] + schur[12, 2] + schur[12, 2, 1] + schur[12, 3] + schur[13, 1, 1] + 2*schur[13, 2] + schur[13, 3] + schur[14, 1] + schur[14, 1, 1] + schur[14, 2] + schur[15, 1] + schur[15, 2] + schur[16, 1] + schur[17, 1] + schur[19],
(6, 5, 4, 3, 2) : schur[4, 4, 2] + schur[5, 3, 1, 1] + schur[5, 3, 3] + schur[5, 4, 1, 1] + schur[6, 1, 1, 1, 1] + schur[6, 3, 1, 1] + schur[6, 3, 2] + schur[6, 4, 1] + schur[6, 4, 2] + schur[6, 5, 1] + schur[6, 6] + schur[7, 2, 1, 1] + schur[7, 3, 1] + schur[7, 3, 1, 1] + schur[7, 3, 2] + 2*schur[7, 4, 1] + schur[7, 5, 1] + schur[7, 6] + schur[8, 1, 1, 1] + schur[8, 2, 1, 1] + schur[8, 2, 2] + 2*schur[8, 3, 1] + schur[8, 3, 2] + schur[8, 4] + 2*schur[8, 4, 1] + schur[8, 5] + schur[8, 6] + schur[9, 1, 1, 1] + schur[9, 2, 1] + schur[9, 2, 1, 1] + 3*schur[9, 3, 1] + schur[9, 4] + schur[9, 4, 1] + schur[9, 5] + schur[10, 1, 1, 1] + 2*schur[10, 2, 1] + schur[10, 2, 2] + schur[10, 3] + 2*schur[10, 3, 1] + 2*schur[10, 4] + schur[10, 5] + schur[11, 1, 1] + schur[11, 1, 1, 1] + 2*schur[11, 2, 1] + 2*schur[11, 3] + schur[11, 3, 1] + schur[11, 4] + schur[12, 1, 1] + schur[12, 2] + 2*schur[12, 2, 1] + 2*schur[12, 3] + schur[12, 4] + 2*schur[13, 1, 1] + schur[13, 2] + schur[13, 2, 1] + schur[13, 3] + schur[14, 1, 1] + 2*schur[14, 2] + schur[14, 3] + schur[15, 1] + schur[15, 1, 1] + schur[15, 2] + schur[16, 1] + schur[16, 2] + schur[17, 1] + schur[18, 1] + schur[20],
(6, 5, 4, 3, 2, 1) : schur[1, 1, 1, 1, 1, 1] + schur[3, 1, 1, 1, 1] + schur[4, 1, 1, 1, 1] + schur[4, 2, 1, 1] + schur[4, 3, 1, 1] + schur[4, 4, 1] + schur[4, 4, 1, 1] + schur[4, 4, 2] + schur[4, 4, 3] + schur[5, 1, 1, 1, 1] + schur[5, 2, 1, 1] + schur[5, 2, 2] + schur[5, 3, 1, 1] + schur[5, 3, 2] + schur[5, 4, 1] + schur[5, 4, 1, 1] + schur[5, 4, 2] + schur[6, 1, 1, 1] + schur[6, 1, 1, 1, 1] + 2*schur[6, 2, 1, 1] + schur[6, 3, 1] + 2*schur[6, 3, 1, 1] + schur[6, 3, 2] + schur[6, 3, 3] + 2*schur[6, 4, 1] + schur[6, 4, 1, 1] + schur[6, 4, 2] + schur[6, 5, 1] + schur[6, 6, 1] + schur[7, 1, 1, 1] + schur[7, 1, 1, 1, 1] + schur[7, 2, 1] + schur[7, 2, 1, 1] + schur[7, 2, 2] + 2*schur[7, 3, 1] + schur[7, 3, 1, 1] + schur[7, 3, 2] + schur[7, 4] + 3*schur[7, 4, 1] + schur[7, 4, 2] + schur[7, 5] + schur[7, 5, 1] + schur[7, 6] + schur[7, 7] + 2*schur[8, 1, 1, 1] + 2*schur[8, 2, 1] + 2*schur[8, 2, 1, 1] + schur[8, 2, 2] + 3*schur[8, 3, 1] + schur[8, 3, 1, 1] + schur[8, 3, 2] + schur[8, 4] + 3*schur[8, 4, 1] + schur[8, 5] + schur[8, 5, 1] + schur[8, 6] + 2*schur[9, 1, 1, 1] + 2*schur[9, 2, 1] + schur[9, 2, 1, 1] + schur[9, 2, 2] + schur[9, 3] + 3*schur[9, 3, 1] + schur[9, 3, 2] + 2*schur[9, 4] + 2*schur[9, 4, 1] + schur[9, 5] + schur[9, 6] + schur[10, 1, 1] + 2*schur[10, 1, 1, 1] + 3*schur[10, 2, 1] + schur[10, 2, 1, 1] + schur[10, 3] + 3*schur[10, 3, 1] + 2*schur[10, 4] + schur[10, 4, 1] + schur[10, 5] + schur[11, 1, 1] + schur[11, 1, 1, 1] + schur[11, 2] + 3*schur[11, 2, 1] + schur[11, 2, 2] + 2*schur[11, 3] + 2*schur[11, 3, 1] + 2*schur[11, 4] + schur[11, 5] + 2*schur[12, 1, 1] + schur[12, 1, 1, 1] + schur[12, 2] + 2*schur[12, 2, 1] + 2*schur[12, 3] + schur[12, 3, 1] + schur[12, 4] + 2*schur[13, 1, 1] + 2*schur[13, 2] + 2*schur[13, 2, 1] + 2*schur[13, 3] + schur[13, 4] + 2*schur[14, 1, 1] + schur[14, 2] + schur[14, 2, 1] + schur[14, 3] + schur[15, 1] + schur[15, 1, 1] + 2*schur[15, 2] + schur[15, 3] + schur[16, 1] + schur[16, 1, 1] + schur[16, 2] + schur[17, 1] + schur[17, 2] + schur[18, 1] + schur[19, 1] + schur[21]
}

## 5,3,1 cover relations

[[[[5, 3, 1], [5, 3, 1]], [[5, 3, 1], [5, 3]]],
 [[[5, 3, 1], [5, 3, 1]], [[5, 3, 1], [5, 2, 1]]],
 [[[5, 3, 1], [5, 3, 1]], [[5, 3, 1], [4, 3, 1]]],
 [[[5, 3, 1], [5, 3]], [[5, 3, 1], [5, 2]]],
 [[[5, 3, 1], [5, 3]], [[5, 3, 1], [4, 3]]],
 [[[5, 3, 1], [5, 2, 1]], [[5, 3, 1], [5, 2]]],
 [[[5, 3, 1], [5, 2, 1]], [[5, 3, 1], [5, 1, 1]]],
 [[[5, 3, 1], [5, 2, 1]], [[5, 3, 1], [4, 2, 1]]],
 [[[5, 3, 1], [5, 2]], [[5, 3, 1], [5, 1]]],
 [[[5, 3, 1], [5, 2]], [[5, 3, 1], [4, 2]]],
 [[[5, 3, 1], [5, 1, 1]], [[5, 3, 1], [5]]],
 [[[5, 3, 1], [5, 1, 1]], [[5, 3, 1], [4, 1, 1]]],
 [[[5, 3, 1], [5, 1]], [[5, 3, 1], [5]]],
 [[[5, 3, 1], [5, 1]], [[5, 3, 1], [4, 1]]],
 [[[5, 3, 1], [5]], [[5, 3, 1], [4]]],
 [[[5, 3, 1], [4, 3, 1]], [[5, 3, 1], [4, 3]]],
 [[[5, 3, 1], [4, 3, 1]], [[5, 3, 1], [4, 2, 1]]],
 [[[5, 3, 1], [4, 3, 1]], [[5, 3, 1], [3, 3, 1]]],
 [[[5, 3, 1], [4, 3]], [[5, 3, 1], [4, 2]]],
 [[[5, 3, 1], [4, 3]], [[5, 3, 1], [3, 3]]],
 [[[5, 3, 1], [4, 2, 1]], [[5, 3, 1], [4, 2]]],
 [[[5, 3, 1], [4, 2, 1]], [[5, 3, 1], [4, 1, 1]]],
 [[[5, 3, 1], [4, 2, 1]], [[5, 3, 1], [3, 2, 1]]],
 [[[5, 3, 1], [4, 2]], [[5, 3, 1], [4, 1]]],
 [[[5, 3, 1], [4, 2]], [[5, 3, 1], [3, 2]]],
 [[[5, 3, 1], [4, 1, 1]], [[5, 3, 1], [4]]],
 [[[5, 3, 1], [4, 1, 1]], [[5, 3, 1], [3, 1, 1]]],
 [[[5, 3, 1], [4, 1]], [[5, 3, 1], [4]]],
 [[[5, 3, 1], [4, 1]], [[5, 3, 1], [3, 1]]],
 [[[5, 3, 1], [4]], [[5, 3, 1], [3]]],
 [[[5, 3, 1], [3, 3, 1]], [[5, 3, 1], [3, 3]]],
 [[[5, 3, 1], [3, 3, 1]], [[5, 3, 1], [2, 2, 1]]],
 [[[5, 3, 1], [3, 3]], [[5, 3, 1], [2, 2]]],
 [[[5, 3, 1], [3, 2, 1]], [[5, 3, 1], [3, 2]]],
 [[[5, 3, 1], [3, 2, 1]], [[5, 3, 1], [3, 1, 1]]],
 [[[5, 3, 1], [3, 2, 1]], [[5, 3, 1], [2, 2, 1]]],
 [[[5, 3, 1], [3, 2]], [[5, 3, 1], [3, 1]]],
 [[[5, 3, 1], [3, 2]], [[5, 3, 1], [2, 2]]],
 [[[5, 3, 1], [3, 1, 1]], [[5, 3, 1], [2, 1, 1]]],
 [[[5, 3, 1], [3, 1, 1]], [[5, 3, 1], [3]]],
 [[[5, 3, 1], [3, 1]], [[5, 3, 1], [3]]],
 [[[5, 3, 1], [3, 1]], [[5, 3, 1], [2, 1]]],
 [[[5, 3, 1], [3]], [[5, 3, 1], [2]]],
 [[[5, 3, 1], [2, 2, 1]], [[5, 3, 1], [2, 2]]],
 [[[5, 3, 1], [2, 2, 1]], [[5, 3, 1], [1, 1, 1]]],
 [[[5, 3, 1], [2, 2]], [[5, 3, 1], [1, 1]]],
 [[[5, 3, 1], [2, 1, 1]], [[5, 3, 1], [1, 1, 1]]],
 [[[5, 3, 1], [2, 1, 1]], [[5, 3, 1], [1]]],
 [[[5, 3, 1], [2, 1]], [[5, 3, 1], [2]]],
 [[[5, 3, 1], [2, 1]], [[5, 3, 1], [1, 1]]],
 [[[5, 3, 1], [2]], [[5, 3, 1], [1]]],
 [[[5, 3, 1], [1, 1, 1]], [[5, 3, 1], []]],
 [[[5, 3, 1], [1, 1]], [[5, 3, 1], []]],
 [[[5, 3, 1], [1]], [[5, 3, 1], []]]]

# cover relations 6,4,2,1

[[[[6, 4, 2, 1], [6, 4, 2, 1]], [[6, 4, 2, 1], [6, 4, 2]]],
 [[[6, 4, 2, 1], [6, 4, 2, 1]], [[6, 4, 2, 1], [6, 4, 1, 1]]],
 [[[6, 4, 2, 1], [6, 4, 2, 1]], [[6, 4, 2, 1], [6, 3, 2, 1]]],
 [[[6, 4, 2, 1], [6, 4, 2, 1]], [[6, 4, 2, 1], [5, 4, 2, 1]]],
 [[[6, 4, 2, 1], [6, 4, 2]], [[6, 4, 2, 1], [6, 4, 1]]],
 [[[6, 4, 2, 1], [6, 4, 2]], [[6, 4, 2, 1], [6, 3, 2]]],
 [[[6, 4, 2, 1], [6, 4, 2]], [[6, 4, 2, 1], [5, 4, 2]]],
 [[[6, 4, 2, 1], [6, 4, 1, 1]], [[6, 4, 2, 1], [6, 4]]],
 [[[6, 4, 2, 1], [6, 4, 1, 1]], [[6, 4, 2, 1], [6, 3, 1, 1]]],
 [[[6, 4, 2, 1], [6, 4, 1, 1]], [[6, 4, 2, 1], [5, 4, 1, 1]]],
 [[[6, 4, 2, 1], [6, 4, 1]], [[6, 4, 2, 1], [6, 4]]],
 [[[6, 4, 2, 1], [6, 4, 1]], [[6, 4, 2, 1], [6, 3, 1]]],
 [[[6, 4, 2, 1], [6, 4, 1]], [[6, 4, 2, 1], [5, 4, 1]]],
 [[[6, 4, 2, 1], [6, 4]], [[6, 4, 2, 1], [6, 3]]],
 [[[6, 4, 2, 1], [6, 4]], [[6, 4, 2, 1], [5, 4]]],
 [[[6, 4, 2, 1], [6, 3, 2, 1]], [[6, 4, 2, 1], [6, 3, 2]]],
 [[[6, 4, 2, 1], [6, 3, 2, 1]], [[6, 4, 2, 1], [6, 3, 1, 1]]],
 [[[6, 4, 2, 1], [6, 3, 2, 1]], [[6, 4, 2, 1], [6, 2, 2, 1]]],
 [[[6, 4, 2, 1], [6, 3, 2, 1]], [[6, 4, 2, 1], [5, 3, 2, 1]]],
 [[[6, 4, 2, 1], [6, 3, 2]], [[6, 4, 2, 1], [6, 3, 1]]],
 [[[6, 4, 2, 1], [6, 3, 2]], [[6, 4, 2, 1], [6, 2, 2]]],
 [[[6, 4, 2, 1], [6, 3, 2]], [[6, 4, 2, 1], [5, 3, 2]]],
 [[[6, 4, 2, 1], [6, 3, 1, 1]], [[6, 4, 2, 1], [6, 3]]],
 [[[6, 4, 2, 1], [6, 3, 1, 1]], [[6, 4, 2, 1], [6, 2, 1, 1]]],
 [[[6, 4, 2, 1], [6, 3, 1, 1]], [[6, 4, 2, 1], [5, 3, 1, 1]]],
 [[[6, 4, 2, 1], [6, 3, 1]], [[6, 4, 2, 1], [6, 3]]],
 [[[6, 4, 2, 1], [6, 3, 1]], [[6, 4, 2, 1], [6, 2, 1]]],
 [[[6, 4, 2, 1], [6, 3, 1]], [[6, 4, 2, 1], [5, 3, 1]]],
 [[[6, 4, 2, 1], [6, 3]], [[6, 4, 2, 1], [6, 2]]],
 [[[6, 4, 2, 1], [6, 3]], [[6, 4, 2, 1], [5, 3]]],
 [[[6, 4, 2, 1], [6, 2, 2, 1]], [[6, 4, 2, 1], [6, 2, 2]]],
 [[[6, 4, 2, 1], [6, 2, 2, 1]], [[6, 4, 2, 1], [6, 1, 1, 1]]],
 [[[6, 4, 2, 1], [6, 2, 2, 1]], [[6, 4, 2, 1], [5, 2, 2, 1]]],
 [[[6, 4, 2, 1], [6, 2, 2]], [[6, 4, 2, 1], [6, 1, 1]]],
 [[[6, 4, 2, 1], [6, 2, 2]], [[6, 4, 2, 1], [5, 2, 2]]],
 [[[6, 4, 2, 1], [6, 2, 1, 1]], [[6, 4, 2, 1], [6, 1, 1, 1]]],
 [[[6, 4, 2, 1], [6, 2, 1, 1]], [[6, 4, 2, 1], [6, 1]]],
 [[[6, 4, 2, 1], [6, 2, 1, 1]], [[6, 4, 2, 1], [5, 2, 1, 1]]],
 [[[6, 4, 2, 1], [6, 2, 1]], [[6, 4, 2, 1], [6, 2]]],
 [[[6, 4, 2, 1], [6, 2, 1]], [[6, 4, 2, 1], [6, 1, 1]]],
 [[[6, 4, 2, 1], [6, 2, 1]], [[6, 4, 2, 1], [5, 2, 1]]],
 [[[6, 4, 2, 1], [6, 2]], [[6, 4, 2, 1], [6, 1]]],
 [[[6, 4, 2, 1], [6, 2]], [[6, 4, 2, 1], [5, 2]]],
 [[[6, 4, 2, 1], [6, 1, 1, 1]], [[6, 4, 2, 1], [6]]],
 [[[6, 4, 2, 1], [6, 1, 1, 1]], [[6, 4, 2, 1], [5, 1, 1, 1]]],
 [[[6, 4, 2, 1], [6, 1, 1]], [[6, 4, 2, 1], [6]]],
 [[[6, 4, 2, 1], [6, 1, 1]], [[6, 4, 2, 1], [5, 1, 1]]],
 [[[6, 4, 2, 1], [6, 1]], [[6, 4, 2, 1], [6]]],
 [[[6, 4, 2, 1], [6, 1]], [[6, 4, 2, 1], [5, 1]]],
 [[[6, 4, 2, 1], [6]], [[6, 4, 2, 1], [5]]],
 [[[6, 4, 2, 1], [5, 4, 2, 1]], [[6, 4, 2, 1], [5, 4, 2]]],
 [[[6, 4, 2, 1], [5, 4, 2, 1]], [[6, 4, 2, 1], [5, 4, 1, 1]]],
 [[[6, 4, 2, 1], [5, 4, 2, 1]], [[6, 4, 2, 1], [5, 3, 2, 1]]],
 [[[6, 4, 2, 1], [5, 4, 2, 1]], [[6, 4, 2, 1], [4, 4, 2, 1]]],
 [[[6, 4, 2, 1], [5, 4, 2]], [[6, 4, 2, 1], [5, 4, 1]]],
 [[[6, 4, 2, 1], [5, 4, 2]], [[6, 4, 2, 1], [5, 3, 2]]],
 [[[6, 4, 2, 1], [5, 4, 2]], [[6, 4, 2, 1], [4, 4, 2]]],
 [[[6, 4, 2, 1], [5, 4, 1, 1]], [[6, 4, 2, 1], [5, 4]]],
 [[[6, 4, 2, 1], [5, 4, 1, 1]], [[6, 4, 2, 1], [5, 3, 1, 1]]],
 [[[6, 4, 2, 1], [5, 4, 1, 1]], [[6, 4, 2, 1], [4, 4, 1, 1]]],
 [[[6, 4, 2, 1], [5, 4, 1]], [[6, 4, 2, 1], [5, 4]]],
 [[[6, 4, 2, 1], [5, 4, 1]], [[6, 4, 2, 1], [5, 3, 1]]],
 [[[6, 4, 2, 1], [5, 4, 1]], [[6, 4, 2, 1], [4, 4, 1]]],
 [[[6, 4, 2, 1], [5, 4]], [[6, 4, 2, 1], [5, 3]]],
 [[[6, 4, 2, 1], [5, 4]], [[6, 4, 2, 1], [4, 4]]],
 [[[6, 4, 2, 1], [5, 3, 2, 1]], [[6, 4, 2, 1], [5, 3, 2]]],
 [[[6, 4, 2, 1], [5, 3, 2, 1]], [[6, 4, 2, 1], [5, 3, 1, 1]]],
 [[[6, 4, 2, 1], [5, 3, 2, 1]], [[6, 4, 2, 1], [5, 2, 2, 1]]],
 [[[6, 4, 2, 1], [5, 3, 2, 1]], [[6, 4, 2, 1], [4, 3, 2, 1]]],
 [[[6, 4, 2, 1], [5, 3, 2]], [[6, 4, 2, 1], [5, 3, 1]]],
 [[[6, 4, 2, 1], [5, 3, 2]], [[6, 4, 2, 1], [5, 2, 2]]],
 [[[6, 4, 2, 1], [5, 3, 2]], [[6, 4, 2, 1], [4, 3, 2]]],
 [[[6, 4, 2, 1], [5, 3, 1, 1]], [[6, 4, 2, 1], [5, 3]]],
 [[[6, 4, 2, 1], [5, 3, 1, 1]], [[6, 4, 2, 1], [5, 2, 1, 1]]],
 [[[6, 4, 2, 1], [5, 3, 1, 1]], [[6, 4, 2, 1], [4, 3, 1, 1]]],
 [[[6, 4, 2, 1], [5, 3, 1]], [[6, 4, 2, 1], [5, 3]]],
 [[[6, 4, 2, 1], [5, 3, 1]], [[6, 4, 2, 1], [5, 2, 1]]],
 [[[6, 4, 2, 1], [5, 3, 1]], [[6, 4, 2, 1], [4, 3, 1]]],
 [[[6, 4, 2, 1], [5, 3]], [[6, 4, 2, 1], [5, 2]]],
 [[[6, 4, 2, 1], [5, 3]], [[6, 4, 2, 1], [4, 3]]],
 [[[6, 4, 2, 1], [5, 2, 2, 1]], [[6, 4, 2, 1], [5, 2, 2]]],
 [[[6, 4, 2, 1], [5, 2, 2, 1]], [[6, 4, 2, 1], [5, 1, 1, 1]]],
 [[[6, 4, 2, 1], [5, 2, 2, 1]], [[6, 4, 2, 1], [4, 2, 2, 1]]],
 [[[6, 4, 2, 1], [5, 2, 2]], [[6, 4, 2, 1], [5, 1, 1]]],
 [[[6, 4, 2, 1], [5, 2, 2]], [[6, 4, 2, 1], [4, 2, 2]]],
 [[[6, 4, 2, 1], [5, 2, 1, 1]], [[6, 4, 2, 1], [5, 1, 1, 1]]],
 [[[6, 4, 2, 1], [5, 2, 1, 1]], [[6, 4, 2, 1], [5, 1]]],
 [[[6, 4, 2, 1], [5, 2, 1, 1]], [[6, 4, 2, 1], [4, 2, 1, 1]]],
 [[[6, 4, 2, 1], [5, 2, 1]], [[6, 4, 2, 1], [5, 2]]],
 [[[6, 4, 2, 1], [5, 2, 1]], [[6, 4, 2, 1], [5, 1, 1]]],
 [[[6, 4, 2, 1], [5, 2, 1]], [[6, 4, 2, 1], [4, 2, 1]]],
 [[[6, 4, 2, 1], [5, 2]], [[6, 4, 2, 1], [5, 1]]],
 [[[6, 4, 2, 1], [5, 2]], [[6, 4, 2, 1], [4, 2]]],
 [[[6, 4, 2, 1], [5, 1, 1, 1]], [[6, 4, 2, 1], [5]]],
 [[[6, 4, 2, 1], [5, 1, 1, 1]], [[6, 4, 2, 1], [4, 1, 1, 1]]],
 [[[6, 4, 2, 1], [5, 1, 1]], [[6, 4, 2, 1], [5]]],
 [[[6, 4, 2, 1], [5, 1, 1]], [[6, 4, 2, 1], [4, 1, 1]]],
 [[[6, 4, 2, 1], [5, 1]], [[6, 4, 2, 1], [5]]],
 [[[6, 4, 2, 1], [5, 1]], [[6, 4, 2, 1], [4, 1]]],
 [[[6, 4, 2, 1], [5]], [[6, 4, 2, 1], [4]]],
 [[[6, 4, 2, 1], [4, 4, 2, 1]], [[6, 4, 2, 1], [4, 4, 2]]],
 [[[6, 4, 2, 1], [4, 4, 2, 1]], [[6, 4, 2, 1], [4, 4, 1, 1]]],
 [[[6, 4, 2, 1], [4, 4, 2, 1]], [[6, 4, 2, 1], [3, 3, 2, 1]]],
 [[[6, 4, 2, 1], [4, 4, 2]], [[6, 4, 2, 1], [4, 4, 1]]],
 [[[6, 4, 2, 1], [4, 4, 2]], [[6, 4, 2, 1], [3, 3, 2]]],
 [[[6, 4, 2, 1], [4, 4, 1, 1]], [[6, 4, 2, 1], [4, 4]]],
 [[[6, 4, 2, 1], [4, 4, 1, 1]], [[6, 4, 2, 1], [3, 3, 1, 1]]],
 [[[6, 4, 2, 1], [4, 4, 1]], [[6, 4, 2, 1], [4, 4]]],
 [[[6, 4, 2, 1], [4, 4, 1]], [[6, 4, 2, 1], [3, 3, 1]]],
 [[[6, 4, 2, 1], [4, 4]], [[6, 4, 2, 1], [3, 3]]],
 [[[6, 4, 2, 1], [4, 3, 2, 1]], [[6, 4, 2, 1], [4, 3, 2]]],
 [[[6, 4, 2, 1], [4, 3, 2, 1]], [[6, 4, 2, 1], [4, 3, 1, 1]]],
 [[[6, 4, 2, 1], [4, 3, 2, 1]], [[6, 4, 2, 1], [4, 2, 2, 1]]],
 [[[6, 4, 2, 1], [4, 3, 2, 1]], [[6, 4, 2, 1], [3, 3, 2, 1]]],
 [[[6, 4, 2, 1], [4, 3, 2]], [[6, 4, 2, 1], [4, 3, 1]]],
 [[[6, 4, 2, 1], [4, 3, 2]], [[6, 4, 2, 1], [4, 2, 2]]],
 [[[6, 4, 2, 1], [4, 3, 2]], [[6, 4, 2, 1], [3, 3, 2]]],
 [[[6, 4, 2, 1], [4, 3, 1, 1]], [[6, 4, 2, 1], [4, 3]]],
 [[[6, 4, 2, 1], [4, 3, 1, 1]], [[6, 4, 2, 1], [4, 2, 1, 1]]],
 [[[6, 4, 2, 1], [4, 3, 1, 1]], [[6, 4, 2, 1], [3, 3, 1, 1]]],
 [[[6, 4, 2, 1], [4, 3, 1]], [[6, 4, 2, 1], [4, 3]]],
 [[[6, 4, 2, 1], [4, 3, 1]], [[6, 4, 2, 1], [4, 2, 1]]],
 [[[6, 4, 2, 1], [4, 3, 1]], [[6, 4, 2, 1], [3, 3, 1]]],
 [[[6, 4, 2, 1], [4, 3]], [[6, 4, 2, 1], [4, 2]]],
 [[[6, 4, 2, 1], [4, 3]], [[6, 4, 2, 1], [3, 3]]],
 [[[6, 4, 2, 1], [4, 2, 2, 1]], [[6, 4, 2, 1], [4, 2, 2]]],
 [[[6, 4, 2, 1], [4, 2, 2, 1]], [[6, 4, 2, 1], [4, 1, 1, 1]]],
 [[[6, 4, 2, 1], [4, 2, 2, 1]], [[6, 4, 2, 1], [3, 2, 2, 1]]],
 [[[6, 4, 2, 1], [4, 2, 2]], [[6, 4, 2, 1], [4, 1, 1]]],
 [[[6, 4, 2, 1], [4, 2, 2]], [[6, 4, 2, 1], [3, 2, 2]]],
 [[[6, 4, 2, 1], [4, 2, 1, 1]], [[6, 4, 2, 1], [4, 1, 1, 1]]],
 [[[6, 4, 2, 1], [4, 2, 1, 1]], [[6, 4, 2, 1], [3, 2, 1, 1]]],
 [[[6, 4, 2, 1], [4, 2, 1, 1]], [[6, 4, 2, 1], [4, 1]]],
 [[[6, 4, 2, 1], [4, 2, 1]], [[6, 4, 2, 1], [4, 2]]],
 [[[6, 4, 2, 1], [4, 2, 1]], [[6, 4, 2, 1], [4, 1, 1]]],
 [[[6, 4, 2, 1], [4, 2, 1]], [[6, 4, 2, 1], [3, 2, 1]]],
 [[[6, 4, 2, 1], [4, 2]], [[6, 4, 2, 1], [4, 1]]],
 [[[6, 4, 2, 1], [4, 2]], [[6, 4, 2, 1], [3, 2]]],
 [[[6, 4, 2, 1], [4, 1, 1, 1]], [[6, 4, 2, 1], [3, 1, 1, 1]]],
 [[[6, 4, 2, 1], [4, 1, 1, 1]], [[6, 4, 2, 1], [4]]],
 [[[6, 4, 2, 1], [4, 1, 1]], [[6, 4, 2, 1], [4]]],
 [[[6, 4, 2, 1], [4, 1, 1]], [[6, 4, 2, 1], [3, 1, 1]]],
 [[[6, 4, 2, 1], [4, 1]], [[6, 4, 2, 1], [4]]],
 [[[6, 4, 2, 1], [4, 1]], [[6, 4, 2, 1], [3, 1]]],
 [[[6, 4, 2, 1], [4]], [[6, 4, 2, 1], [3]]],
 [[[6, 4, 2, 1], [3, 3, 2, 1]], [[6, 4, 2, 1], [3, 3, 2]]],
 [[[6, 4, 2, 1], [3, 3, 2, 1]], [[6, 4, 2, 1], [3, 3, 1, 1]]],
 [[[6, 4, 2, 1], [3, 3, 2, 1]], [[6, 4, 2, 1], [2, 2, 2, 1]]],
 [[[6, 4, 2, 1], [3, 3, 2]], [[6, 4, 2, 1], [3, 3, 1]]],
 [[[6, 4, 2, 1], [3, 3, 2]], [[6, 4, 2, 1], [2, 2, 2]]],
 [[[6, 4, 2, 1], [3, 3, 1, 1]], [[6, 4, 2, 1], [3, 3]]],
 [[[6, 4, 2, 1], [3, 3, 1, 1]], [[6, 4, 2, 1], [2, 2, 1, 1]]],
 [[[6, 4, 2, 1], [3, 3, 1]], [[6, 4, 2, 1], [3, 3]]],
 [[[6, 4, 2, 1], [3, 3, 1]], [[6, 4, 2, 1], [2, 2, 1]]],
 [[[6, 4, 2, 1], [3, 3]], [[6, 4, 2, 1], [2, 2]]],
 [[[6, 4, 2, 1], [3, 2, 2, 1]], [[6, 4, 2, 1], [3, 2, 2]]],
 [[[6, 4, 2, 1], [3, 2, 2, 1]], [[6, 4, 2, 1], [2, 2, 2, 1]]],
 [[[6, 4, 2, 1], [3, 2, 2, 1]], [[6, 4, 2, 1], [2, 1, 1, 1]]],
 [[[6, 4, 2, 1], [3, 2, 2]], [[6, 4, 2, 1], [2, 2, 2]]],
 [[[6, 4, 2, 1], [3, 2, 2]], [[6, 4, 2, 1], [2, 1, 1]]],
 [[[6, 4, 2, 1], [3, 2, 1, 1]], [[6, 4, 2, 1], [3, 1, 1, 1]]],
 [[[6, 4, 2, 1], [3, 2, 1, 1]], [[6, 4, 2, 1], [2, 2, 1, 1]]],
 [[[6, 4, 2, 1], [3, 2, 1, 1]], [[6, 4, 2, 1], [2, 1]]],
 [[[6, 4, 2, 1], [3, 2, 1]], [[6, 4, 2, 1], [3, 2]]],
 [[[6, 4, 2, 1], [3, 2, 1]], [[6, 4, 2, 1], [3, 1, 1]]],
 [[[6, 4, 2, 1], [3, 2, 1]], [[6, 4, 2, 1], [2, 2, 1]]],
 [[[6, 4, 2, 1], [3, 2]], [[6, 4, 2, 1], [3, 1]]],
 [[[6, 4, 2, 1], [3, 2]], [[6, 4, 2, 1], [2, 2]]],
 [[[6, 4, 2, 1], [3, 1, 1, 1]], [[6, 4, 2, 1], [2, 1, 1, 1]]],
 [[[6, 4, 2, 1], [3, 1, 1, 1]], [[6, 4, 2, 1], [2]]],
 [[[6, 4, 2, 1], [3, 1, 1]], [[6, 4, 2, 1], [3]]],
 [[[6, 4, 2, 1], [3, 1, 1]], [[6, 4, 2, 1], [2, 1, 1]]],
 [[[6, 4, 2, 1], [3, 1]], [[6, 4, 2, 1], [3]]],
 [[[6, 4, 2, 1], [3, 1]], [[6, 4, 2, 1], [2, 1]]],
 [[[6, 4, 2, 1], [3]], [[6, 4, 2, 1], [2]]],
 [[[6, 4, 2, 1], [2, 2, 2, 1]], [[6, 4, 2, 1], [2, 2, 2]]],
 [[[6, 4, 2, 1], [2, 2, 2, 1]], [[6, 4, 2, 1], [1, 1, 1, 1]]],
 [[[6, 4, 2, 1], [2, 2, 2]], [[6, 4, 2, 1], [1, 1, 1]]],
 [[[6, 4, 2, 1], [2, 2, 1, 1]], [[6, 4, 2, 1], [1, 1, 1, 1]]],
 [[[6, 4, 2, 1], [2, 2, 1, 1]], [[6, 4, 2, 1], [1, 1]]],
 [[[6, 4, 2, 1], [2, 2, 1]], [[6, 4, 2, 1], [2, 2]]],
 [[[6, 4, 2, 1], [2, 2, 1]], [[6, 4, 2, 1], [1, 1, 1]]],
 [[[6, 4, 2, 1], [2, 2]], [[6, 4, 2, 1], [1, 1]]],
 [[[6, 4, 2, 1], [2, 1, 1, 1]], [[6, 4, 2, 1], [1, 1, 1, 1]]],
 [[[6, 4, 2, 1], [2, 1, 1, 1]], [[6, 4, 2, 1], [1]]],
 [[[6, 4, 2, 1], [2, 1, 1]], [[6, 4, 2, 1], [1, 1, 1]]],
 [[[6, 4, 2, 1], [2, 1, 1]], [[6, 4, 2, 1], [1]]],
 [[[6, 4, 2, 1], [2, 1]], [[6, 4, 2, 1], [2]]],
 [[[6, 4, 2, 1], [2, 1]], [[6, 4, 2, 1], [1, 1]]],
 [[[6, 4, 2, 1], [2]], [[6, 4, 2, 1], [1]]],
 [[[6, 4, 2, 1], [1, 1, 1, 1]], [[6, 4, 2, 1], []]],
 [[[6, 4, 2, 1], [1, 1, 1]], [[6, 4, 2, 1], []]],
 [[[6, 4, 2, 1], [1, 1]], [[6, 4, 2, 1], []]],
 [[[6, 4, 2, 1], [1]], [[6, 4, 2, 1], []]]]


# cover relations 3,2,2,1

[[[[3, 2, 2, 1], [3, 2, 2, 1]], [[3, 2, 2, 1], [3, 2, 2]]],
 [[[3, 2, 2, 1], [3, 2, 2, 1]], [[3, 2, 2, 1], [3, 2, 1, 1]]],
 [[[3, 2, 2, 1], [3, 2, 2, 1]], [[3, 2, 2, 1], [2, 2, 2, 1]]],
 [[[3, 2, 2, 1], [3, 2, 2]], [[3, 2, 2, 1], [3, 2]]],
 [[[3, 2, 2, 1], [3, 2, 2]], [[3, 2, 2, 1], [2, 2, 2]]],
 [[[3, 2, 2, 1], [3, 2, 1, 1]], [[3, 2, 2, 1], [3, 2, 1]]],
 [[[3, 2, 2, 1], [3, 2, 1, 1]], [[3, 2, 2, 1], [3, 1, 1, 1]]],
 [[[3, 2, 2, 1], [3, 2, 1, 1]], [[3, 2, 2, 1], [2, 2, 1, 1]]],
 [[[3, 2, 2, 1], [3, 2, 1]], [[3, 2, 2, 1], [3, 2]]],
 [[[3, 2, 2, 1], [3, 2, 1]], [[3, 2, 2, 1], [3, 1, 1]]],
 [[[3, 2, 2, 1], [3, 2, 1]], [[3, 2, 2, 1], [2, 2, 1]]],
 [[[3, 2, 2, 1], [3, 2]], [[3, 2, 2, 1], [3]]],
 [[[3, 2, 2, 1], [3, 2]], [[3, 2, 2, 1], [2, 1]]],
 [[[3, 2, 2, 1], [3, 1, 1, 1]], [[3, 2, 2, 1], [3, 1, 1]]],
 [[[3, 2, 2, 1], [3, 1, 1, 1]], [[3, 2, 2, 1], [1, 1, 1, 1]]],
 [[[3, 2, 2, 1], [3, 1, 1]], [[3, 2, 2, 1], [3, 1]]],
 [[[3, 2, 2, 1], [3, 1, 1]], [[3, 2, 2, 1], [1, 1, 1]]],
 [[[3, 2, 2, 1], [3, 1]], [[3, 2, 2, 1], [3]]],
 [[[3, 2, 2, 1], [3, 1]], [[3, 2, 2, 1], [1]]],
 [[[3, 2, 2, 1], [3]], [[3, 2, 2, 1], []]],
 [[[3, 2, 2, 1], [2, 2, 2, 1]], [[3, 2, 2, 1], [2, 2, 2]]],
 [[[3, 2, 2, 1], [2, 2, 2, 1]], [[3, 2, 2, 1], [2, 1, 1, 1]]],
 [[[3, 2, 2, 1], [2, 2, 2]], [[3, 2, 2, 1], [2, 1]]],
 [[[3, 2, 2, 1], [2, 2, 1, 1]], [[3, 2, 2, 1], [2, 2, 1]]],
 [[[3, 2, 2, 1], [2, 2, 1, 1]], [[3, 2, 2, 1], [2, 1, 1, 1]]],
 [[[3, 2, 2, 1], [2, 2, 1]], [[3, 2, 2, 1], [2, 2]]],
 [[[3, 2, 2, 1], [2, 2, 1]], [[3, 2, 2, 1], [2, 1, 1]]],
 [[[3, 2, 2, 1], [2, 2]], [[3, 2, 2, 1], [2]]],
 [[[3, 2, 2, 1], [2, 1, 1, 1]], [[3, 2, 2, 1], [2, 1, 1]]],
 [[[3, 2, 2, 1], [2, 1, 1, 1]], [[3, 2, 2, 1], [1, 1, 1, 1]]],
 [[[3, 2, 2, 1], [2, 1, 1]], [[3, 2, 2, 1], [2, 1]]],
 [[[3, 2, 2, 1], [2, 1, 1]], [[3, 2, 2, 1], [1, 1, 1]]],
 [[[3, 2, 2, 1], [2, 1]], [[3, 2, 2, 1], [2]]],
 [[[3, 2, 2, 1], [2, 1]], [[3, 2, 2, 1], [1, 1]]],
 [[[3, 2, 2, 1], [2]], [[3, 2, 2, 1], []]],
 [[[3, 2, 2, 1], [1, 1, 1, 1]], [[3, 2, 2, 1], [1, 1, 1]]],
 [[[3, 2, 2, 1], [1, 1, 1]], [[3, 2, 2, 1], [1]]],
 [[[3, 2, 2, 1], [1, 1]], [[3, 2, 2, 1], [1]]],
 [[[3, 2, 2, 1], [1]], [[3, 2, 2, 1], []]]]


# cover relations of 2,1,1

[[[[2, 1, 1], [2, 1, 1]], [[2, 1, 1], [2, 1]]],
 [[[2, 1, 1], [2, 1, 1]], [[2, 1, 1], [1, 1, 1]]],
 [[[2, 1, 1], [2, 1]], [[2, 1, 1], [2]]],
 [[[2, 1, 1], [2, 1]], [[2, 1, 1], [1, 1]]],
 [[[2, 1, 1], [2]], [[2, 1, 1], []]],
 [[[2, 1, 1], [1, 1, 1]], [[2, 1, 1], [1]]],
 [[[2, 1, 1], [1, 1]], [[2, 1, 1], [1]]],
 [[[2, 1, 1], [1]], [[2, 1, 1], []]]]

# cover relations of 3,2,1

[[[[3, 2, 1], [3, 2, 1]], [[3, 2, 1], [3, 2]]],
 [[[3, 2, 1], [3, 2, 1]], [[3, 2, 1], [3, 1, 1]]],
 [[[3, 2, 1], [3, 2, 1]], [[3, 2, 1], [2, 2, 1]]],
 [[[3, 2, 1], [3, 2]], [[3, 2, 1], [3, 1]]],
 [[[3, 2, 1], [3, 2]], [[3, 2, 1], [2, 2]]],
 [[[3, 2, 1], [3, 1, 1]], [[3, 2, 1], [3]]],
 [[[3, 2, 1], [3, 1, 1]], [[3, 2, 1], [2, 1, 1]]],
 [[[3, 2, 1], [3, 1]], [[3, 2, 1], [3]]],
 [[[3, 2, 1], [3, 1]], [[3, 2, 1], [2, 1]]],
 [[[3, 2, 1], [3]], [[3, 2, 1], [1]]],
 [[[3, 2, 1], [2, 2, 1]], [[3, 2, 1], [2, 2]]],
 [[[3, 2, 1], [2, 2, 1]], [[3, 2, 1], [1, 1, 1]]],
 [[[3, 2, 1], [2, 2]], [[3, 2, 1], [1, 1]]],
 [[[3, 2, 1], [2, 1, 1]], [[3, 2, 1], [1, 1, 1]]],
 [[[3, 2, 1], [2, 1, 1]], [[3, 2, 1], [1]]],
 [[[3, 2, 1], [2, 1]], [[3, 2, 1], [2]]],
 [[[3, 2, 1], [2, 1]], [[3, 2, 1], [1, 1]]],
 [[[3, 2, 1], [2]], [[3, 2, 1], [1]]],
 [[[3, 2, 1], [1, 1, 1]], [[3, 2, 1], []]],
 [[[3, 2, 1], [1, 1]], [[3, 2, 1], []]],
 [[[3, 2, 1], [1]], [[3, 2, 1], []]]]

# cover relations of 2,1,1,1

[[[[2, 1, 1, 1], [2, 1, 1, 1]], [[2, 1, 1, 1], [2, 1, 1]]],
 [[[2, 1, 1, 1], [2, 1, 1, 1]], [[2, 1, 1, 1], [1, 1, 1, 1]]],
 [[[2, 1, 1, 1], [2, 1, 1]], [[2, 1, 1, 1], [2, 1]]],
 [[[2, 1, 1, 1], [2, 1, 1]], [[2, 1, 1, 1], [1, 1, 1]]],
 [[[2, 1, 1, 1], [2, 1]], [[2, 1, 1, 1], [2]]],
 [[[2, 1, 1, 1], [2, 1]], [[2, 1, 1, 1], [1]]],
 [[[2, 1, 1, 1], [2]], [[2, 1, 1, 1], []]],
 [[[2, 1, 1, 1], [1, 1, 1, 1]], [[2, 1, 1, 1], [1, 1]]],
 [[[2, 1, 1, 1], [1, 1, 1]], [[2, 1, 1, 1], [1, 1]]],
 [[[2, 1, 1, 1], [1, 1]], [[2, 1, 1, 1], [1]]],
 [[[2, 1, 1, 1], [1]], [[2, 1, 1, 1], []]]]

# cover relations 2,1,1,1,1

[[[[2, 1, 1, 1, 1], [2, 1, 1, 1, 1]], [[2, 1, 1, 1, 1], [2, 1, 1, 1]]],
 [[[2, 1, 1, 1, 1], [2, 1, 1, 1, 1]], [[2, 1, 1, 1, 1], [1, 1, 1, 1, 1]]],
 [[[2, 1, 1, 1, 1], [2, 1, 1, 1]], [[2, 1, 1, 1, 1], [2, 1, 1]]],
 [[[2, 1, 1, 1, 1], [2, 1, 1, 1]], [[2, 1, 1, 1, 1], [1, 1, 1, 1]]],
 [[[2, 1, 1, 1, 1], [2, 1, 1]], [[2, 1, 1, 1, 1], [2, 1]]],
 [[[2, 1, 1, 1, 1], [2, 1, 1]], [[2, 1, 1, 1, 1], [1, 1]]],
 [[[2, 1, 1, 1, 1], [2, 1]], [[2, 1, 1, 1, 1], [2]]],
 [[[2, 1, 1, 1, 1], [2, 1]], [[2, 1, 1, 1, 1], [1]]],
 [[[2, 1, 1, 1, 1], [2]], [[2, 1, 1, 1, 1], []]],
 [[[2, 1, 1, 1, 1], [1, 1, 1, 1, 1]], [[2, 1, 1, 1, 1], [1, 1, 1]]],
 [[[2, 1, 1, 1, 1], [1, 1, 1, 1]], [[2, 1, 1, 1, 1], [1, 1, 1]]],
 [[[2, 1, 1, 1, 1], [1, 1, 1]], [[2, 1, 1, 1, 1], [1, 1]]],
 [[[2, 1, 1, 1, 1], [1, 1]], [[2, 1, 1, 1, 1], [1]]],
 [[[2, 1, 1, 1, 1], [1]], [[2, 1, 1, 1, 1], []]]]

# cover relations of 2,2,1,1,1

[[[[2, 2, 1, 1, 1], [2, 2, 1, 1, 1]], [[2, 2, 1, 1, 1], [2, 2, 1, 1]]],
 [[[2, 2, 1, 1, 1], [2, 2, 1, 1, 1]], [[2, 2, 1, 1, 1], [2, 1, 1, 1, 1]]],
 [[[2, 2, 1, 1, 1], [2, 2, 1, 1]], [[2, 2, 1, 1, 1], [2, 2, 1]]],
 [[[2, 2, 1, 1, 1], [2, 2, 1, 1]], [[2, 2, 1, 1, 1], [2, 1, 1, 1]]],
 [[[2, 2, 1, 1, 1], [2, 2, 1]], [[2, 2, 1, 1, 1], [2, 2]]],
 [[[2, 2, 1, 1, 1], [2, 2, 1]], [[2, 2, 1, 1, 1], [2, 1]]],
 [[[2, 2, 1, 1, 1], [2, 2]], [[2, 2, 1, 1, 1], [2]]],
 [[[2, 2, 1, 1, 1], [2, 1, 1, 1, 1]], [[2, 2, 1, 1, 1], [2, 1, 1]]],
 [[[2, 2, 1, 1, 1], [2, 1, 1, 1, 1]], [[2, 2, 1, 1, 1], [1, 1, 1, 1, 1]]],
 [[[2, 2, 1, 1, 1], [2, 1, 1, 1]], [[2, 2, 1, 1, 1], [2, 1, 1]]],
 [[[2, 2, 1, 1, 1], [2, 1, 1, 1]], [[2, 2, 1, 1, 1], [1, 1, 1, 1]]],
 [[[2, 2, 1, 1, 1], [2, 1, 1]], [[2, 2, 1, 1, 1], [2, 1]]],
 [[[2, 2, 1, 1, 1], [2, 1, 1]], [[2, 2, 1, 1, 1], [1, 1, 1]]],
 [[[2, 2, 1, 1, 1], [2, 1]], [[2, 2, 1, 1, 1], [2]]],
 [[[2, 2, 1, 1, 1], [2, 1]], [[2, 2, 1, 1, 1], [1]]],
 [[[2, 2, 1, 1, 1], [2]], [[2, 2, 1, 1, 1], []]],
 [[[2, 2, 1, 1, 1], [1, 1, 1, 1, 1]], [[2, 2, 1, 1, 1], [1, 1]]],
 [[[2, 2, 1, 1, 1], [1, 1, 1, 1]], [[2, 2, 1, 1, 1], [1, 1, 1]]],
 [[[2, 2, 1, 1, 1], [1, 1, 1]], [[2, 2, 1, 1, 1], [1, 1]]],
 [[[2, 2, 1, 1, 1], [1, 1]], [[2, 2, 1, 1, 1], [1]]],
 [[[2, 2, 1, 1, 1], [1]], [[2, 2, 1, 1, 1], []]]]


# cover relations of 3,2,1,1

[[[[3, 2, 1, 1], [3, 2, 1, 1]], [[3, 2, 1, 1], [3, 2, 1]]],
 [[[3, 2, 1, 1], [3, 2, 1, 1]], [[3, 2, 1, 1], [3, 1, 1, 1]]],
 [[[3, 2, 1, 1], [3, 2, 1, 1]], [[3, 2, 1, 1], [2, 2, 1, 1]]],
 [[[3, 2, 1, 1], [3, 2, 1]], [[3, 2, 1, 1], [3, 2]]],
 [[[3, 2, 1, 1], [3, 2, 1]], [[3, 2, 1, 1], [3, 1, 1]]],
 [[[3, 2, 1, 1], [3, 2, 1]], [[3, 2, 1, 1], [2, 2, 1]]],
 [[[3, 2, 1, 1], [3, 2]], [[3, 2, 1, 1], [3]]],
 [[[3, 2, 1, 1], [3, 2]], [[3, 2, 1, 1], [2, 2]]],
 [[[3, 2, 1, 1], [3, 1, 1, 1]], [[3, 2, 1, 1], [3, 1]]], #
 [[[3, 2, 1, 1], [3, 1, 1, 1]], [[3, 2, 1, 1], [1, 1, 1, 1]]],
 [[[3, 2, 1, 1], [3, 1, 1]], [[3, 2, 1, 1], [3, 1]]],
 [[[3, 2, 1, 1], [3, 1, 1]], [[3, 2, 1, 1], [1, 1, 1]]],
 [[[3, 2, 1, 1], [3, 1]], [[3, 2, 1, 1], [3]]],
 [[[3, 2, 1, 1], [3, 1]], [[3, 2, 1, 1], [1]]],
 [[[3, 2, 1, 1], [3]], [[3, 2, 1, 1], []]],
 [[[3, 2, 1, 1], [2, 2, 1, 1]], [[3, 2, 1, 1], [2, 1, 1]]], #
 [[[3, 2, 1, 1], [2, 2, 1, 1]], [[3, 2, 1, 1], [2, 1, 1, 1]]],
 [[[3, 2, 1, 1], [2, 2, 1]], [[3, 2, 1, 1], [2, 2]]],
 [[[3, 2, 1, 1], [2, 2, 1]], [[3, 2, 1, 1], [2, 1, 1]]],
 [[[3, 2, 1, 1], [2, 2]], [[3, 2, 1, 1], [2]]],
 [[[3, 2, 1, 1], [2, 1, 1, 1]], [[3, 2, 1, 1], [2,1]]], #
 [[[3, 2, 1, 1], [2, 1, 1, 1]], [[3, 2, 1, 1], [1, 1, 1, 1]]],
 [[[3, 2, 1, 1], [2, 1, 1]], [[3, 2, 1, 1], [2, 1]]],
 [[[3, 2, 1, 1], [2, 1, 1]], [[3, 2, 1, 1], [1, 1, 1]]],
 [[[3, 2, 1, 1], [2, 1]], [[3, 2, 1, 1], [2]]],
 [[[3, 2, 1, 1], [2, 1]], [[3, 2, 1, 1], [1, 1]]],
 [[[3, 2, 1, 1], [2]], [[3, 2, 1, 1], []]],
 [[[3, 2, 1, 1], [1, 1, 1, 1]], [[3, 2, 1, 1], [1]]], #
 [[[3, 2, 1, 1], [1, 1, 1]], [[3, 2, 1, 1], [1]]], #
 [[[3, 2, 1, 1], [1, 1]], [[3, 2, 1, 1], [1]]],
 [[[3, 2, 1, 1], [1]], [[3, 2, 1, 1], []]]]


# cover relations of 4,2,1 LATTICE GOOD DISTRIBUTION
# Works also with the schur / e expression

[[[[4, 2, 1], [4, 2, 1]], [[4, 2, 1], [4, 2]]],
 [[[4, 2, 1], [4, 2, 1]], [[4, 2, 1], [4, 1, 1]]],
 [[[4, 2, 1], [4, 2, 1]], [[4, 2, 1], [3, 2, 1]]],
 [[[4, 2, 1], [4, 2]], [[4, 2, 1], [4, 1]]],
 [[[4, 2, 1], [4, 2]], [[4, 2, 1], [3, 1]]],
 [[[4, 2, 1], [4, 1, 1]], [[4, 2, 1], [4]]],
 [[[4, 2, 1], [4, 1, 1]], [[4, 2, 1], [2, 1, 1]]],
 [[[4, 2, 1], [4, 1]], [[4, 2, 1], [4]]],
 [[[4, 2, 1], [4, 1]], [[4, 2, 1], [2, 1]]],
 [[[4, 2, 1], [4]], [[4, 2, 1], [1]]],
 [[[4, 2, 1], [3, 2, 1]], [[4, 2, 1], [3, 2]]],
 [[[4, 2, 1], [3, 2, 1]], [[4, 2, 1], [3, 1, 1]]],
 [[[4, 2, 1], [3, 2, 1]], [[4, 2, 1], [2, 2, 1]]],
 [[[4, 2, 1], [3, 2]], [[4, 2, 1], [3, 1]]],
 [[[4, 2, 1], [3, 2]], [[4, 2, 1], [2, 2]]],
 [[[4, 2, 1], [3, 1, 1]], [[4, 2, 1], [3]]],
 [[[4, 2, 1], [3, 1, 1]], [[4, 2, 1], [2, 1, 1]]],
 [[[4, 2, 1], [3, 1]], [[4, 2, 1], [3]]],
 [[[4, 2, 1], [3, 1]], [[4, 2, 1], [2, 1]]],
 [[[4, 2, 1], [3]], [[4, 2, 1], [2]]],
 [[[4, 2, 1], [2, 2, 1]], [[4, 2, 1], [2, 2]]],
 [[[4, 2, 1], [2, 2, 1]], [[4, 2, 1], [1, 1, 1]]],
 [[[4, 2, 1], [2, 2]], [[4, 2, 1], [1, 1]]],
 [[[4, 2, 1], [2, 1, 1]], [[4, 2, 1], [1, 1, 1]]],
 [[[4, 2, 1], [2, 1, 1]], [[4, 2, 1], [1]]],
 [[[4, 2, 1], [2, 1]], [[4, 2, 1], [2]]],
 [[[4, 2, 1], [2, 1]], [[4, 2, 1], [1, 1]]],
 [[[4, 2, 1], [2]], [[4, 2, 1], [1]]],
 [[[4, 2, 1], [1, 1, 1]], [[4, 2, 1], []]],
 [[[4, 2, 1], [1, 1]], [[4, 2, 1], []]],
 [[[4, 2, 1], [1]], [[4, 2, 1], []]]]

# cover relations of 4,3,1

[[[[4, 3, 1], [4, 3, 1]], [[4, 3, 1], [4, 3]]],
 [[[4, 3, 1], [4, 3, 1]], [[4, 3, 1], [4, 2, 1]]],
 [[[4, 3, 1], [4, 3, 1]], [[4, 3, 1], [3, 3, 1]]],
 [[[4, 3, 1], [4, 3]], [[4, 3, 1], [4, 1]]],
 [[[4, 3, 1], [4, 3]], [[4, 3, 1], [3, 3]]],
 [[[4, 3, 1], [4, 2, 1]], [[4, 3, 1], [4, 2]]],
 [[[4, 3, 1], [4, 2, 1]], [[4, 3, 1], [4, 1, 1]]],
 [[[4, 3, 1], [4, 2, 1]], [[4, 3, 1], [3, 2, 1]]],
 [[[4, 3, 1], [4, 2]], [[4, 3, 1], [4, 1]]],
 [[[4, 3, 1], [4, 1, 1]], [[4, 3, 1], [4]]],
 [[[4, 3, 1], [4, 1]], [[4, 3, 1], [4]]],
 [[[4, 3, 1], [3, 3, 1]], [[4, 3, 1], [3, 3]]],
 [[[4, 3, 1], [3, 2, 1]], [[4, 3, 1], [3, 2]]],
 [[[4, 3, 1], [3, 2, 1]], [[4, 3, 1], [3, 1, 1]]],
 [[[4, 3, 1], [3, 2, 1]], [[4, 3, 1], [2, 2, 1]]],
 [[[4, 3, 1], [3, 2]], [[4, 3, 1], [3, 1]]],
 [[[4, 3, 1], [3, 2]], [[4, 3, 1], [2, 2]]],
 [[[4, 3, 1], [3, 1, 1]], [[4, 3, 1], [3]]],
 [[[4, 3, 1], [3, 1, 1]], [[4, 3, 1], [2, 1, 1]]],
 [[[4, 3, 1], [3, 1]], [[4, 3, 1], [3]]],
 [[[4, 3, 1], [3, 1]], [[4, 3, 1], [2, 1]]],
 [[[4, 3, 1], [3]], [[4, 3, 1], [2]]],
 [[[4, 3, 1], [2, 2, 1]], [[4, 3, 1], [2, 2]]],
 [[[4, 3, 1], [2, 2, 1]], [[4, 3, 1], [1, 1, 1]]],
 [[[4, 3, 1], [2, 2]], [[4, 3, 1], [1, 1]]],
 [[[4, 3, 1], [2, 1, 1]], [[4, 3, 1], [1, 1, 1]]],
 [[[4, 3, 1], [2, 1, 1]], [[4, 3, 1], [1]]],
 [[[4, 3, 1], [2, 1]], [[4, 3, 1], [2]]],
 [[[4, 3, 1], [2, 1]], [[4, 3, 1], [1, 1]]],
 [[[4, 3, 1], [2]], [[4, 3, 1], [1]]],
 [[[4, 3, 1], [1, 1, 1]], [[4, 3, 1], []]],
 [[[4, 3, 1], [1, 1]], [[4, 3, 1], []]],
 [[[4, 3, 1], [1]], [[4, 3, 1], []]]]


# Not a lattice but good distribution

[[[[4, 3, 1], [4, 3, 1]], [[4, 3, 1], [4, 3]]],
 [[[4, 3, 1], [4, 3, 1]], [[4, 3, 1], [4, 2, 1]]],
 [[[4, 3, 1], [4, 3, 1]], [[4, 3, 1], [3, 3, 1]]],
 [[[4, 3, 1], [4, 3]], [[4, 3, 1], [4, 1]]],
 [[[4, 3, 1], [4, 3]], [[4, 3, 1], [3, 3]]],
 [[[4, 3, 1], [4, 2, 1]], [[4, 3, 1], [4, 2]]],
 [[[4, 3, 1], [4, 2, 1]], [[4, 3, 1], [4, 1, 1]]],
 [[[4, 3, 1], [4, 2, 1]], [[4, 3, 1], [3, 2, 1]]],
 [[[4, 3, 1], [4, 2]], [[4, 3, 1], [4, 1]]],
 [[[4, 3, 1], [4, 2]], [[4, 3, 1], [3, 2]]],
 [[[4, 3, 1], [4, 1]], [[4, 3, 1], [4]]],
 [[[4, 3, 1], [4, 1]], [[4, 3, 1], [2, 1]]],
 [[[4, 3, 1], [3, 3, 1]], [[4, 3, 1], [4,1,1]]], # also works with 3,1 ?
 [[[4, 3, 1], [3, 3, 1]], [[4, 3, 1], [3, 3]]],
 [[[4, 3, 1], [4, 1, 1]], [[4, 3, 1], [4]]],
 [[[4, 3, 1], [4, 1, 1]], [[4, 3, 1], [2, 1]]],
 [[[4, 3, 1], [4]], [[4, 3, 1], [2]]],
 [[[4, 3, 1], [3, 3]], [[4, 3, 1], [1, 1]]],
 [[[4, 3, 1], [3, 2, 1]], [[4, 3, 1], [3, 2]]],
 [[[4, 3, 1], [3, 2, 1]], [[4, 3, 1], [3, 1, 1]]],
 [[[4, 3, 1], [3, 2, 1]], [[4, 3, 1], [2, 2, 1]]],
 [[[4, 3, 1], [3, 2]], [[4, 3, 1], [3, 1]]],
 [[[4, 3, 1], [3, 2]], [[4, 3, 1], [2, 2]]],
 [[[4, 3, 1], [3, 1, 1]], [[4, 3, 1], [3]]],
 [[[4, 3, 1], [3, 1, 1]], [[4, 3, 1], [2, 1, 1]]],
 [[[4, 3, 1], [3, 1]], [[4, 3, 1], [3]]],
 [[[4, 3, 1], [3, 1]], [[4, 3, 1], [2, 1]]],
 [[[4, 3, 1], [3]], [[4, 3, 1], [2]]],
 [[[4, 3, 1], [2, 2, 1]], [[4, 3, 1], [2, 2]]],
 [[[4, 3, 1], [2, 2, 1]], [[4, 3, 1], [1, 1, 1]]],
 [[[4, 3, 1], [2, 2]], [[4, 3, 1], [1, 1]]],
 [[[4, 3, 1], [2, 1, 1]], [[4, 3, 1], [1, 1, 1]]],
 [[[4, 3, 1], [2, 1, 1]], [[4, 3, 1], [1]]],
 [[[4, 3, 1], [2, 1]], [[4, 3, 1], [2]]],
 [[[4, 3, 1], [2, 1]], [[4, 3, 1], [1, 1]]],
 [[[4, 3, 1], [2]], [[4, 3, 1], [1]]],
 [[[4, 3, 1], [1, 1, 1]], [[4, 3, 1], []]],
 [[[4, 3, 1], [1, 1]], [[4, 3, 1], []]],
 [[[4, 3, 1], [1]], [[4, 3, 1], []]]]

