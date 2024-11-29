
from sage.arith.misc import gcd
from sage.arith.misc import multinomial
from sage.categories.enumerated_sets import EnumeratedSets
from sage.categories.finite_enumerated_sets import FiniteEnumeratedSets
from sage.combinat.composition import Composition
from sage.combinat.partition import Partition
from sage.combinat.partition import Partitions
from sage.combinat.posets.lattices import LatticePoset
from sage.combinat.skew_partition import SkewPartition
from sage.combinat.sf.sf import SymmetricFunctions
from sage.combinat.tableau import Tableau
from sage.functions.other import ceil
from sage.misc.cachefunc import cached_method
from sage.misc.cachefunc import cached_function
from sage.misc.inherit_comparison import InheritComparisonClasscallMetaclass
from sage.misc.latex import latex
from sage.misc.lazy_attribute import lazy_attribute
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.rational_field import QQ
from sage.sets.disjoint_union_enumerated_sets import DisjointUnionEnumeratedSets
from sage.sets.family import Family
from sage.sets.non_negative_integers import NonNegativeIntegers
from sage.structure.element import Element
from sage.structure.parent import Parent
from sage.structure.richcmp import richcmp
from sage.structure.unique_representation import UniqueRepresentation
from sage.symbolic.ring import SR

Partitions.options(latex="young_diagram", convention="french")

class TriangularPartition(Element,
        metaclass=InheritComparisonClasscallMetaclass):

    @staticmethod
    def __classcall_private__(cls, *args, **opts):
        r"""
        Ensure that triangular partition created by the enumerated sets and
        directly are the same and that they are instances of
        :class:`TriangularPartition`.

        TESTS::

            sage: tp = TriangularPartition([3,1])
            sage: tp.parent()
            Triangular Partitions
            sage: type(tp)
            <class '__main__.TriangularPartitions_all_with_category.element_class'>
            sage: tp2 = TriangularPartitions()([3,1])
            sage: tp2.parent() is tp.parent()
            True
            sage: type(tp2) is type(tp)
            True
        """
        P = TriangularPartitions_all()
        return P.element_class(P, *args, **opts)

    def __init__(self, parent, partition):
        r"""
        TESTS::

            sage: TriangularPartition([3,1]).parent()
            Triangular Partitions
        """
        Element.__init__(self, parent)
        self._partition = Partition(partition)

    def _repr_(self):
        r"""
        TESTS::

            sage: TriangularPartition([3,1])
            [3, 1]
        """
        return str(self._partition)

    def _richcmp_(self, other, op):
        if not isinstance(other, TriangularPartition):
            return NotImplemented
        return richcmp(self._partition, other._partition,op)

    def __hash__(self):
        return hash(self._partition)

    def __iter__(self):
        r"""
        TESTS::

            sage: list(TriangularPartition([3,1]))
            [3, 1]
        """
        return iter(self._partition)

    def __getitem__(self, key):
        r"""
        TESTS::

            sage: TriangularPartition([3,1])[0]
            3
        """
        return self._partition.__getitem__(key)

    def _latex_(self):
        r"""
        TESTS::

            sage: s = latex(TriangularPartition([3,1]))
        """
        #return latex(self.partition())
        latex.add_package_to_preamble_if_available("tikz")
        return self._grid_latex_()

    def _grid_latex_(self):
        maxx = max(self)
        maxy = self.length()
        x,y = 0,0
        path = f"({x},{y})"
        for p in self:
            path+= f"--({p},{y})"
            y+=1
            path+= f"--({p},{y})"
        path+=f"--(0,{y})--cycle;\n"
        s = "\\begin{tikzpicture}\n"
        s+= f"\\draw[ultra thick] " + path
        s+= f"\\draw[step=1cm,gray,very thin] (0,0) grid ({maxx},{maxy});\n"
        s+= "\\end{tikzpicture}"
        return s

    def partition(self):
        r"""
        EXAMPLES::

            sage: TriangularPartition([3,1]).partition()
            [3, 1]
            sage: TriangularPartition([3,1]).partition().parent()
            Partitions
        """
        return self._partition

    def size(self):
        r"""
        EXAMPLES::

            sage: TriangularPartition([3,1]).size()
            4
        """
        return self.partition().size()

    def length(self):
        r"""
        EXAMPLES::

            sage: TriangularPartition([3,1]).length()
            2
        """
        return self.partition().length()

    def sub_partitions(self):
        r"""
        EXAMPLES::

            sage: list(TriangularPartition([3,2]).sub_partitions())
            [[3, 2], [3, 1], [2, 2], [3], [2, 1], [1, 1], [2], [1], []]
        """
        l0 = set([self])
        l1 = set()
        while len(l0) > 0:
            for el in l0:
                yield el
                l1.update(TriangularPartition(p) for p in el.partition().down())
            l0 = l1
            l1 = set()

    def conjugate(self):
        r"""
        EXAMPLES::

            sage: TriangularPartition([3,1]).conjugate()
            [2, 1, 1]
        """
        return TriangularPartition(self.partition().conjugate())


    def cells(self):
        r"""
        EXAMPLES::

            sage: TriangularPartition([3,1]).cells()
            [(0, 0), (0, 1), (0, 2), (1, 0)]
        """
        return self._partition.cells()

    def cell_min_slope(self,i,j):
        r"""
        EXAMPLES::

            sage: TriangularPartition([3,1]).cell_min_slope(0,0)
            1/4
        """
        p = self._partition
        return p.leg_length(i,j) / p.hook_length(i,j)

    def cell_max_slope(self,i,j):
        r"""
        EXAMPLES::

            sage: TriangularPartition([3,1]).cell_max_slope(0,0)
            1/2
        """
        p = self._partition
        return (p.leg_length(i,j)+1) / p.hook_length(i,j)

    def cell_mean_slope(self,i,j):
        r"""
        EXAMPLES::

            sage: TriangularPartition([3,1]).cell_mean_slope(0,0)
            3/8
        """
        return (self.cell_min_slope(i,j) + self.cell_max_slope(i,j))/2

    @cached_method
    def global_min_slope(self):
        r"""
        EXAMPLES::

            sage: TriangularPartition([3,1]).global_min_slope()
            1/4
        """
        return max((self.cell_min_slope(i,j) for i,j in self._partition.cells()), default = 0)

    @cached_method
    def global_max_slope(self):
        r"""
        EXAMPLES::

            sage: TriangularPartition([3,1]).global_max_slope()
            1/2
        """
        return min((self.cell_max_slope(i,j) for i,j in self._partition.cells()), default = 1)

    @cached_method
    def mean_slope(self):
        r"""
        EXAMPLES::

            sage: TriangularPartition([3,1]).mean_slope()
            3/8
        """
        return (self.global_min_slope() + self.global_max_slope())/2

    def is_triangular(self):
        r"""
        EXAMPLES::

            sage: TriangularPartition([3,1]).is_triangular()
            True
            sage: TriangularPartition([2,2]).is_triangular()
            False
        """
        return self.global_min_slope() < self.global_max_slope()

    def min_slope_cells(self):
        r"""
        EXAMPLES::

            sage: list(TriangularPartition([3,1]).min_slope_cells())
            [(0, 0)]
        """
        for i,j in self.cells():
            if self.cell_min_slope(i,j) == self.global_min_slope():
                yield (i,j)

    def max_slope_cells(self):
        r"""
        EXAMPLES::

            sage: list(TriangularPartition([3,1]).max_slope_cells())
            [(0, 0), (0, 1)]
        """
        for i,j in self.cells():
            if self.cell_max_slope(i,j) == self.global_max_slope():
                yield (i,j)

    def min_slopes(self):
        r"""
        EXAMPLES::

            sage: TriangularPartition([3,1]).min_slopes()
            {(0, 0): 1/4, (0, 1): 0, (0, 2): 0, (1, 0): 0}
        """
        return {(i,j): self.cell_min_slope(i,j) for i,j in self.cells()}

    def max_slopes(self):
        r"""
        EXAMPLES::

            sage: TriangularPartition([3,1]).max_slopes()
            {(0, 0): 1/2, (0, 1): 1/2, (0, 2): 1, (1, 0): 1}
        """
        return {(i,j): self.cell_max_slope(i,j) for i,j in self.cells()}

    def up(self):
        r"""
        EXAMPLES::

            sage: list(TriangularPartition([3,1]).up())
            [[4, 1], [3, 2]]
        """
        for p in self._partition.up():
            pp = TriangularPartition(p)
            if pp.is_triangular():
                yield pp

    def down(self):
        r"""
        EXAMPLES::

            sage: list(TriangularPartition([3,1]).down())
            [[2, 1], [3]]
        """
        for p in self._partition.down():
            pp = TriangularPartition(p)
            if pp.is_triangular():
                yield pp

    def corners(self):
        r"""
        EXAMPLES::

            sage: TriangularPartition([3,1]).corners()
            [(0, 2), (1, 0)]
        """
        return self._partition.corners()

    def horizontal_rectangle_cells(self,c1, c2):
        return [c for c in self.cells() if c[0] > c2[0] and c[0] <= c1[0] and c[1] > c1[1] and c[1] < c2[1]]

    def vertical_rectangle_cells(self,c1,c2):
        return [c for c in self.cells() if c[0] > c2[0] and c[0] < c1[0] and c[1] > c1[1] and c[1] <= c2[1]]

    def triangular_corners(self):
        r"""
        EXAMPLES::

            sage: list(TriangularPartition([3,1]).triangular_corners())
            [(0, 2), (1, 0)]
        """
        for i,j in self.corners():
            p = TriangularPartition(self.partition().remove_cell(i))
            if p.is_triangular():
                yield i,j

    def is_sim_cell(self,i,j,tau):
        r"""
        EXAMPLES::

            sage: lam = TriangularPartition([1,1])
            sage: tau = TriangularPartition([3,1])
            sage: lam.is_sim_cell(0,0,tau)
            False
        """
        return self.cell_min_slope(i,j) < tau.mean_slope() and tau.mean_slope() <= self.cell_max_slope(i,j)

    def sim_cells(self,tau):
        r"""
        EXAMPLES::

            sage: lam = TriangularPartition([1,1])
            sage: tau = TriangularPartition([3,1])
            sage: list(lam.sim_cells(tau))
            [(1, 0)]
        """
        for i,j in self.cells():
            if self.is_sim_cell(i,j,tau):
                yield (i,j)

    def is_similar(self,tau):
        r"""
        EXAMPLES::

            sage: lam = TriangularPartition([1,1])
            sage: tau = TriangularPartition([3,1])
            sage: lam.is_similar(tau)
            False
        """
        return all(self.is_sim_cell(i,j,tau) for i,j in self.cells())

    def contains(self, other):
        r"""
        EXAMPLES::

            sage: lam = TriangularPartition([1,1])
            sage: tau = TriangularPartition([3,1])
            sage: tau.contains(lam)
            True
        """
        return self._partition.contains(other._partition)

    @cached_method
    def interior(self):
        r"""
        EXAMPLES::

            sage: TriangularPartition([3,1]).interior()
            [2]
        """
        return TriangularPartitions.young_meet(self.down())

    @cached_method
    def exterior(self):
        r"""
        EXAMPLES::

            sage: TriangularPartition([3,1]).exterior()
            [4, 2]
        """
        return TriangularPartitions.young_join(self.up())

    @cached_method
    def diagonal(self):
        r"""
        EXAMPLES::

            sage: TriangularPartition([3,1]).diagonal()
            [3, 1] / [2]
        """
        return SkewPartition([self.partition(),self.interior().partition()])


    def decompose(self, cell):
        r"""
        EXAMPLES::

            sage: TriangularPartition([4,3,1]).decompose((1,2))
            ([1], [1], (1, 2))
        """
        y,x = cell
        return TriangularPartition(self[y+1:]), TriangularPartition([v - x - 1 for v in self[:y]]), cell

    def up_diagonal(self,line):
        r"""
        EXAMPLES::

            sage: TriangularPartition([4,3,1]).up_diagonal(1)
            [1] / [1]
        """
        d = self.diagonal()
        return SkewPartition([d[0][line+1:],d[1][line+1:]])

    def diagonal_decompose(self, cell):
        r"""
        EXAMPLES::

            sage: TriangularPartition([4,3,1]).diagonal_decompose((1,2))
            ([1] / [1], [1], (1, 2))
        """
        y,x = cell
        return self.up_diagonal(y), TriangularPartition([v - x - 1 for v in self[:y]]), cell

    def compose(self, cell, other):
        r"""
        EXAMPLES::

            sage: t = TriangularPartition([1])
            sage: t.compose((1,2),t)
            [4, 3, 1]
        """
        y, x = cell
        p2 = list(other._partition)
        p2.extend([0] * (y - len(p2)))
        return TriangularPartition([v + x + 1 for v in p2] + [x+1] + list(self))

    def decompositions(self):
        r"""
        EXAMPLES::

            sage: list(TriangularPartition([4,3,1]).decompositions())
            [([1], [1], (1, 2))]
        """
        for cell in self.diagonal().cells():
            yield self.decompose(cell)

    def diagonal_decompositions(self):
        r"""
        EXAMPLES::

            sage: list(TriangularPartition([4,3,1]).diagonal_decompositions())
            [([1] / [1], [1], (1, 2))]
        """
        for cell in self.diagonal().cells():
            yield self.diagonal_decompose(cell)

    def qAreaEnumerator(self, q = None):
        r"""
        EXAMPLES::

            sage: TriangularPartition([4,3,1]).qAreaEnumerator()
            q^8 + q^7 + 2*q^6 + 3*q^5 + 4*q^4 + 4*q^3 + 4*q^2 + 3*q + 1
        """
        if q is None:
            K = PolynomialRing(QQ,"q")
            q = K.gen()
        return qAreaEnumerator(self, q)


    def Aqt(self, gens = None, tableau = None):
        r"""
        EXAMPLES::

            sage: TriangularPartition([4,3,1]).Aqt()
            q^8 + q^7*t + q^6*t^2 + q^5*t^3 + q^4*t^4 + q^3*t^5 + q^2*t^6 + q*t^7 + t^8 + q^6*t + q^5*t^2 + q^4*t^3 + q^3*t^4 + q^2*t^5 + q*t^6 + q^5*t + 2*q^4*t^2 + 2*q^3*t^3 + 2*q^2*t^4 + q*t^5
        """
        if gens is None:
            K = PolynomialRing(QQ,"q,t")
            gens = K.gens()
        q,t = gens
        return sum(q**tdp.area() * t**tdp.sim(tableau) for tdp in self.triangular_dyck_paths())

    def Aqt_schur2parts(self):
        r"""
        EXAMPLES::

            sage: tp = TriangularPartition([3,2])
            sage: tp.Aqt_schur2parts()
            s[3, 1] + s[5]
        """
        pol = self.Aqt()
        q,t = pol.parent().gens()
        Sym = SymmetricFunctions(QQ)
        Schur = Sym.Schur()
        return sum(coeff * Schur(p) for p,coeff in list(Schur.from_polynomial(pol)) if p.length() <= 2)

    def parking_symmetric_enumeration(self):
        r"""
        EXAMPLES::

            sage: TriangularPartition([4,3,1]).parking_symmetric_enumeration()
            (q^2+q+1)*e[1, 1, 1, 1] + (q^5+2*q^4+4*q^3+2*q^2+2*q)*e[2, 1, 1] + (q^6+q^4+q^2)*e[2, 2] + (q^7+q^6+2*q^5+q^4)*e[3, 1] + q^8*e[4]
        """
        Symq = SymmetricFunctions(PolynomialRing(QQ,"q"))
        q = Symq.base_ring().gens()[0]
        e = Symq.elementary()
        return parkingFunctionSymmetricEnumerator(self, self.length()+1, q,e)

    def number_of_parking_functions(self):
        r"""
        EXAMPLES::

            sage: TriangularPartition([4,3,1]).number_of_parking_functions()
            243
        """
        f = self.parking_symmetric_enumeration()
        return sum(c(1) * multinomial(list(p)) for p,c in f)

    def triangular_dyck_paths(self):
        r"""
        EXAMPLES::

            sage: TriangularPartition([4,3,1]).triangular_dyck_paths()
            Triangular Dyck Paths of [4, 3, 1]
        """
        return TriangularDyckPaths(self)

    @cached_method
    def path_tamari_lattice(self):
        r"""
        EXAMPLES::

            sage: TriangularPartition([4,3,1]).path_tamari_lattice()
            Finite lattice containing 23 elements
        """
        return LatticePoset({k:list(k.path_tamari_up()) for k in self.triangular_dyck_paths()})

    def triangular_tamari_poset(self):
        return Poset({k:list(k.triangular_tamari_up()) for k in self.triangular_dyck_paths()})

    def slope_tamari_poset(self):
        return Poset({k:list(k.slope_tamari_up()) for k in self.triangular_dyck_paths()})


    def immediate_deficit_tamari_poset(self):
        return Poset({k:list(k.immediate_deficit_tamari_up()) for k in self.triangular_dyck_paths()})

    def triangular_schur_poset(self):
        return Poset({k:list(k.triangular_schur_up()) for k in self.triangular_dyck_paths()})

    def area_sim_distribution(self, tableau = None):
        r"""
        EXAMPLES::

            sage: TriangularPartition([3,1]).area_sim_distribution() == {(4, 0): 1, (3, 1): 1, (2, 2): 1, (1, 2): 1, (2, 1): 1, (1, 3): 1, (0, 4): 1}
            True
        """
        d = {}
        for dp in self.triangular_dyck_paths():
            q,t = dp.area(), dp.sim(tableau)
            d[(q,t)] = d.get((q,t),0)+1
        return d

    def top_down_area_sim_distribution(self):
        return self.area_sim_distribution(tableau = self.top_down_standard_tableau())

    def path_distance_sim_distribution_intervals(self):
        r"""
        EXAMPLES::

            sage: TriangularPartition([2]).path_distance_sim_distribution_intervals() == {(0, 2): 1, (1, 1): 1, (2, 0): 1, (0, 1): 1, (1, 0): 1, (0, 0): 1}
            True
        """
        L = self.path_tamari_lattice()
        d = {}
        for dp1, dp2 in L.relations():
            distance = dp1.path_distance(dp2)
            sim = dp2.sim()
            d[(distance, sim)] = d.get((distance, sim),0) + 1
        return d

    def triangular_distance_sim_distribution_intervals(self):
        L = self.triangular_tamari_poset()
        d = {}
        for dp1, dp2 in L.relations():
            distance = dp1.triangular_distance(dp2)
            sim = dp2.sim()
            d[(distance, sim)] = d.get((distance, sim),0) + 1
        return d

    def triangular_min_slope(self):
        t = self.global_min_slope()
        return -t / (1 - t)

    def triangular_max_slope(self):
        t = self.global_max_slope()
        return -t / (1 - t)

    def triangular_min_line(self):
        for i,j in self.cells():
            if self.cell_min_slope(i,j) == self.global_min_slope():
                x = self[i] + 1
                y = i+1
                a = self.triangular_min_slope()
                b = y - a*x
                return a, b

    def triangular_max_line(self):
        for i,j in self.cells():
            if self.cell_max_slope(i,j) == self.global_max_slope():
                x = self[i]
                y = i+1
                a = self.triangular_max_slope()
                b = y - a*x
                return a, b

    def is_max_triangle_path(self,a,b):
        f = lambda x: -a/b * x + a
        return all(i+1 <= f(self[i]) and i+1 > f(self[i]+1) for i in range(self.length())) and self.length() + 1 > f(1)

    def rectangular_slopes(self):
        r"""
        EXAMPLES::

            sage: list(TriangularPartition([3,1]).rectangular_slopes())
            [(3, 5)]
        """
        if self.length() == 1:
            amin = 1
            amax = 3
            bmin = self[0] + 1
            bmax = 3*self[0]
        elif self.conjugate().length() == 1:
            amin = self.length() + 1
            amax = 3*self.length()
            bmin = 1
            bmax = 3
        else:
            lmin, cmin = self.triangular_min_line()
            lmax, cmax = self.triangular_max_line()
            amin = ceil(cmin)
            amax = ceil(cmax)
            bmin = ceil(-cmax/lmax)
            bmax = ceil(-cmin/lmin)
        for a in range(amin,amax):
            for b in range(bmin,bmax):
                if self.is_max_triangle_path(a,b):
                    yield a,b

    def rational_slopes(self):
        r"""
        EXAMPLES::

            sage: list(TriangularPartition([3,1]).rational_slopes())
            [(3, 5)]
        """
        for a,b in self.rectangular_slopes():
            if gcd(a,b) == 1:
                yield a,b

    def is_rectangular(self):
        r"""
        EXAMPLES::

            sage: TriangularPartition([3,1]).is_rectangular()
            True
        """
        for a,b in self.rectangular_slopes():
            return True
        return False

    def is_rational(self):
        r"""
        EXAMPLES::

            sage: TriangularPartition([3,1]).is_rational()
            True
        """
        for a,b in self.rational_slopes():
            return True
        return False

    # @cached_method
    # def natural_corner(self):
        # # if self.size() == 0:
            # # return None
        # # corners = list(self.triangular_corners())
        # # if self.mean_slope() >= 1/2:
            # # return min(corners, key = lambda c:c[1])
        # # else:
            # # return min(corners, key = lambda c:c[0])
        # mu = self.natural_down()
        # if mu is None:
            # return None
        # return SkewPartition([self.partition(),mu.partition()]).cells()[0]

    # @cached_method
    # def natural_down(self):
        # if self.size() == 0:
            # return None
        # #return TriangularPartition(self.partition().remove_cell(self.natural_corner()[0]))
        # for mu in self.down():
            # if mu.is_similar(self):
                # return mu

    @cached_method
    def similar_down(self, mu):
        r"""
        EXAMPLES::

            sage: tau = TriangularPartition([3,1])
            sage: tau.similar_down(tau)
            [2, 1]
            sage: TriangularPartition([2,1]).similar_down(tau)
            [2]
        """
        if self.size() == 0:
            return None
        for d in self.down():
            if d.is_similar(mu):
                return d

    @cached_method
    def similar_corner(self,mu):
        r"""
        EXAMPLES::

            sage: tau = TriangularPartition([3,1])
            sage: tau.similar_corner(tau)
            (0, 2)
        """
        d = self.similar_down(mu)
        if d is None:
            return None
        return SkewPartition([self.partition(),d.partition()]).cells()[0]

    @cached_method
    def triangular_standard_labels(self):
        r"""
        EXAMPLES::

            sage: TriangularPartition([3,1]).triangular_standard_labels() == {(0, 2): 4, (1, 0): 3, (0, 1): 2, (0, 0): 1}
            True
        """
        labels = {}
        p = self
        while p.size() > 0:
            d = p.similar_down(self)
            labels[p.similar_corner(self)] = p.size()
            p = d
        return labels

        # mu = self.natural_down()
        # labels = dict(mu.triangular_standard_labels())
        # labels[self.natural_corner()] = self.size()
        # return labels

    @cached_method
    def top_down_standard_labels(self):
        r"""
        EXAMPLES::

            sage: TriangularPartition([3,1]).top_down_standard_labels() == {(1, 0): 4, (0, 2): 3, (0, 1): 2, (0, 0): 1}
            True
        """
        labels = {}
        l = list(self)
        v = self.size()
        while v > 0:
            for i in range(len(l)-1,-1,-1):
                l[i] -= 1
                labels[(i,l[i])] = v
                v -=1
                if l[i] == 0:
                    l.pop()

        return labels

    @cached_method
    def top_down_standard_tableau(self):
        r"""
        EXAMPLES::

            sage: TriangularPartition([3,1]).top_down_standard_tableau()
            [[1, 2, 3], [4]]
        """
        t = [[0] * k for k in self]
        l = self.top_down_standard_labels()
        for i,j in self.cells():
            t[i][j] = l[(i,j)]
        return Tableau(t)

    @cached_method
    def triangular_standard_tableau(self):
        r"""
        EXAMPLES::

            sage: TriangularPartition([3,1]).triangular_standard_tableau()
            [[1, 2, 4], [3]]
        """
        t = [[0] * k for k in self]
        l = self.triangular_standard_labels()
        for i,j in self.cells():
            t[i][j] = l[(i,j)]
        return Tableau(t)

    def standard_tableau(self, labels):
        t = [[0] * k for k in self]
        l = labels
        for i,j in self.cells():
            t[i][j] = l[(i,j)]
        return Tableau(t)

    def is_sim_sym(self, tableau):
        d_triangular = self.area_sim_distribution()
        return self.area_sim_distribution(tableau = tableau) == d

    def sim_sym_tableaux(self):
        d_triangular = self.area_sim_distribution()
        for t in self.partition().standard_tableaux():
            d = self.area_sim_distribution(t)
            if d == d_triangular:
                yield t

    def diagonal_orientation(self):
        L = self.triangular_standard_labels()
        s = set(L.values())
        p = self
        orientation = None
        while p.size() > 0:
            diagonal = p.diagonal()

            ld = [L[c] for c in diagonal.cells()]
            s.difference_update(ld)
            if not all(y > x for y in ld for x in s):
                return None
            if diagonal.size() > 1:
                dc = diagonal.cells()
                dc.sort(key = lambda c: c[1])
                if orientation is None:
                    if all(L[dc[i]] < L[dc[i+1]] for i in range(len(dc)-1)):
                        orientation = 1
                    elif all(L[dc[i]] > L[dc[i+1]] for i in range(len(dc)-1)):
                        orientation = -1
                    else:
                        return None
                else:
                    if orientation == 1 and any(L[dc[i]] > L[dc[i+1]] for i in range(len(dc)-1)):
                        return None
                    if orientation == -1 and any(L[dc[i]] < L[dc[i+1]] for i in range(len(dc)-1)):
                        return None

            p = p.interior()
        return orientation

    def is_diagonal_oriented(self):
        L = self.triangular_standard_labels()
        s = set(L.values())
        p = self
        while p.size() > 0:
            diagonal = p.diagonal()

            ld = [L[c] for c in diagonal.cells()]
            s.difference_update(ld)
            if not all(y > x for y in ld for x in s):
                return False
            dc = diagonal.cells()
            dc.sort(key = lambda c: c[1])
            if any(L[dc[i]] < L[dc[i+1]] for i in range(len(dc)-1)) and any(L[dc[i]] > L[dc[i+1]] for i in range(len(dc)-1)):
                return False
            p = p.interior()
        return True

    # def is_corner_extremal(self):
        # if self.size() == 0:
            # return True
        # L = self.triangular_standard_labels()
        # max_corner =max(self.corners(), key = lambda x: L[x])
        # C = sorted(self.corners(), key = lambda x: x[0])
        # return max_corner == C[0] or max_corner == C[-1]

    def is_unambiguous(self):
        L = self.triangular_standard_labels()
        conj = self.conjugate()
        for c in self.cells():
            if c[0] > 0 and c[1] > 0:
                if any(L[(i,min(c[1]-1,self[i]-1))] < L[c] for i in range(c[0]+1,self.length())) and any(L[(min(c[0]-1,conj[j]-1),j)] < L[c] for j in range(c[1]+1,conj.length())):
                    return False
        return True

    def slope_enumerator(self, step = 0.01):
        t = self.global_min_slope() + step
        while t < self.global_max_slope():
            yield t
            t+= step

    def slope_label(self, slope):
        t = slope
        values = {}
        distincs = set()
        for a,b in self.cells():
            x = b+1
            y = a+1
            v = x*t + y * (1-t)
            values[(a,b)] =v
            distincs.add(v)
        if len(distincs) == len(values):
            l = list(values)
            l.sort(key = lambda x: values[x])
            labels = {c:i+1 for i,c in enumerate(l)}
            return labels

    def slope_labels(self, step = 0.01):
        for t in self.slope_enumerator(step):
            labels = self.slope_label(t)
            if labels:
                yield labels

    def slope_labels_distincts(self, step = 0.01):
        r"""
        EXAMPLES::

            sage: list(TriangularPartition([3,1]).slope_labels_distincts()) == [{(0, 0): 1, (0, 1): 2, (0, 2): 3, (1, 0): 4}, {(0, 0): 1, (0, 1): 2, (1, 0): 3, (0, 2): 4}]
            True
        """
        tableaux = set()
        for t in self.slope_enumerator(step):
            labels = self.slope_label(t)
            if labels:
                tab = self.standard_tableau(labels)
                if not tab in tableaux:
                    tableaux.add(tab)
                    yield labels


    def slope_tableaux(self, step = 0.01):
        r"""
        EXAMPLES::

            sage: TriangularPartition([3,1]).slope_tableaux()
            {[[1, 2, 3], [4]], [[1, 2, 4], [3]]}
        """
        return set(self.standard_tableau(labels) for labels in self.slope_labels(step))

        return tableaux

    def slope_tableau(self, slope):
        labels = self.slope_label(slope)
        if labels:
            return self.standard_tableau(labels)


    # @cached_method
    # def diagonal_labels(self):
        # if self.size() == 0:
            # return {}
        # labels = {}
        # d = self.diagonal()
        # labels_int = self.interior().diagonal_labels()
        # for k in labels_int:
            # labels[k] = labels_int[k]+1
        # for cell in d.cells():
            # labels[cell] = 1
        # return labels

    # def diagonal_label(self,cell):
        # return self.diagonal_labels()[cell]

    # @cached_method
    # def diagonal_standard_labels(self):
        # if self.size() == 0:
            # return {}
        # labels = {}
        # d = self.diagonal()
        # labels_int = self.interior().diagonal_standard_labels()
        # for k in labels_int:
            # labels[k] = labels_int[k]+d.size()
        # C = d.cells()
        # if self.mean_slope() >= 1/2:
            # C.sort(key = lambda c:c[1])
        # else:
            # C.sort(key = lambda c:c[0])
        # i = 1
        # for cell in C:
            # labels[cell] = i
            # i+=1
        # return labels

    # def diagonal_standard_label(self,cell):
        # return self.diagonal_standard_labels()[cell]

    # def diagonal_tableau(self):
        # t = [[0] * k for k in self]
        # l = self.diagonal_labels()
        # for i,j in self.cells():
            # t[i][j] = l[(i,j)]
        # return Tableau(t)

    # def diagonal_standard_tableau(self):
        # t = [[0] * k for k in self]
        # l = self.diagonal_standard_labels()
        # for i,j in self.cells():
            # t[i][j] = l[(i,j)]
        # return Tableau(t)


@cached_function
def qAreaEnumerator(partition, q):
    if partition.size() == 0:
        return 1
    return q**partition.diagonal().size() * partition.interior().qAreaEnumerator(q) + sum( q**skewalpha.size()* TriangularPartition(skewalpha[1]).qAreaEnumerator(q) * beta.qAreaEnumerator(q) for skewalpha,beta, cell in partition.diagonal_decompositions())

@cached_function
def parkingFunctionSymmetricEnumerator(partition, n, q, e):
    if partition.size() == 0:
        return e[n]
    return q**partition.diagonal().size() * parkingFunctionSymmetricEnumerator(partition.interior(),n,q,e) + sum( q**skewalpha.size() * parkingFunctionSymmetricEnumerator(TriangularPartition(skewalpha[1]), n - cell[0] - 1, q, e) * parkingFunctionSymmetricEnumerator(beta,cell[0]+1,q,e) for skewalpha,beta, cell in partition.diagonal_decompositions())



class TriangularPartitions(UniqueRepresentation, Parent):

    def __call__(self, *args, **keywords):
        r"""
        TESTS::

            sage: TriangularPartitions()([4,3,1])
            [4, 3, 1]
        """

        if isinstance(args[0], TriangularPartition):
            return args[0]

        return super(TriangularPartitions, self).__call__(*args, **keywords)

    @staticmethod
    def __classcall_private__(cls, n=None, **keywords):
        r"""
        TESTS::

            sage: TriangularPartitions()
            Triangular Partitions
            sage: TriangularPartitions(5)
            Triangular Partitions of the integer 5
        """

        if n is None:
            return TriangularPartitions_all()

        return TriangularPartitions_size(n, **keywords)

    @staticmethod
    def young_meet(l):
        r"""
        EXAMPLES::

            sage: p1 = TriangularPartition([2])
            sage: p2 = TriangularPartition([1,1])
            sage: TriangularPartitions().young_meet([p1,p2])
            [1]
        """
        l = list(l)
        if len(l) == 0:
            return TriangularPartition([])
        l0 = {l[0]}
        l1 = set()
        while True:
            for tp in l0:
                if all(tp2.contains(tp) for tp2 in l):
                    return tp
                l1.update(tp.down())
            l0 = l1
            l1 = set()

    @staticmethod
    def young_join(l):
        r"""
        EXAMPLES::

            sage: p1 = TriangularPartition([2])
            sage: p2 = TriangularPartition([1,1])
            sage: TriangularPartitions().young_join([p1,p2])
            [2, 1]
        """
        l = list(l)
        if len(l) == 0:
            return TriangularPartition([])
        l0 = {l[0]}
        l1 = set()
        while True:
            for tp in l0:
                if all(tp.contains(tp2) for tp2 in l):
                    return tp
                l1.update(tp.up())
            l0 = l1
            l1 = set()



class TriangularPartitions_all(DisjointUnionEnumeratedSets, TriangularPartitions):

    def __init__(self):
        DisjointUnionEnumeratedSets.__init__(
            self, Family(NonNegativeIntegers(), TriangularPartitions_size),
            facade=True, keepkey=False,
            category=EnumeratedSets())

    def _repr_(self):
        r"""
        TESTS::

            sage: TriangularPartitions()
            Triangular Partitions
        """
        return "Triangular Partitions"

    def _element_constructor_(self, partition):
        r"""
        EXAMPLES::

            sage: TP = TriangularPartitions()
            sage: TP([4,3,1])
            [4, 3, 1]

        """
        return self.element_class(self, partition)

    def __contains__(self, x):
        r"""
        TESTS::

            sage: TP = TriangularPartitions()
            sage: 1 in TP
            False
            sage: TP([4,3,1]) in TP
            True

        """
        return isinstance(x, self.element_class)

    Element = TriangularPartition

class TriangularPartitions_size(TriangularPartitions):
    r"""
    The enumerated set of Trianglar partitions of a given size
    """
    def __init__(self, size, **keywords):
        # sage: for i in range(15): TestSuite(TriangularPartitions(i)).run()
        r"""
        TESTS::

            sage: TP5 = TriangularPartitions(5)
            sage: assert TP5 is TriangularPartitions(5)
        """
        super(TriangularPartitions_size, self).__init__(category=FiniteEnumeratedSets())

        self._size = size
        self._partitions = Partitions(size, **keywords)

    def _repr_(self):
        r"""
        TESTS::

            sage: TriangularPartitions(5)
            Triangular Partitions of the integer 5
        """
        return "Triangular {}".format(self._partitions)

    def __contains__(self, x):
        r"""
        TESTS::

            sage: TP5 = TriangularPartitions(5)
            sage: 1 in TP5
            False
            sage: TP5([4,1]) in TP5
            True
        """
        return isinstance(x, self.element_class) and x.size() == self._size

    @lazy_attribute
    def _parent_for(self):
        r"""
        The parent of the element generated by ``self``.

        TESTS::

            sage: TP5 = TriangularPartitions(5)
            sage: TP5._parent_for
            Triangular Partitions
        """
        return TriangularPartitions_all()

    # This is needed because this is a facade parent via DisjointUnionEnumeratedSets
    @lazy_attribute
    def element_class(self):
        r"""
        TESTS::

            sage: TP5 = TriangularPartitions(5)
            sage: TP5.element_class
            <class '__main__.TriangularPartitions_all_with_category.element_class'>
            sage: TP5.first().__class__ == TriangularPartitions().first().__class__
            True
        """
        return self._parent_for.element_class

    def _element_constructor_(self, partition):
        r"""
        EXAMPLES::

            sage: TP5 = TriangularPartitions(5)
            sage: TP5([4,1])
            [4, 1]
        """
        return self.element_class(self, partition)

    def __iter__(self):
        r"""
        TESTS::

            sage: TP5 = TriangularPartitions(5)
            sage: list(TP5)
            [[5], [4, 1], [3, 2], [2, 2, 1], [2, 1, 1, 1], [1, 1, 1, 1, 1]]
        """
        for p in self._partitions:
            pp = TriangularPartition(p)
            if pp.is_triangular():
                yield pp

    def rectangular_partitions(self):
        r"""
        EXAMPLES::

            sage: TP5 = TriangularPartitions(5)
            sage: list(TP5.rectangular_partitions())
            [[5], [1, 1, 1, 1, 1]]
            sage: TP6 = TriangularPartitions(6)
            sage: list(TP6.rectangular_partitions())
            [[6], [4, 2], [3, 2, 1], [2, 2, 1, 1], [1, 1, 1, 1, 1, 1]]
        """
        for tr in self:
            if tr.is_rectangular():
                yield tr

    def rational_partitions(self):
        r"""
        EXAMPLES::

            sage: TP5 = TriangularPartitions(5)
            sage: list(TP5.rational_partitions())
            [[5], [1, 1, 1, 1, 1]]
            sage: TP6 = TriangularPartitions(6)
            sage: list(TP6.rational_partitions())
            [[6], [4, 2], [3, 2, 1], [2, 2, 1, 1], [1, 1, 1, 1, 1, 1]]
        """
        for tr in self:
            if tr.is_rational():
                yield tr

class TriangularDyckPath(Element,
        metaclass=InheritComparisonClasscallMetaclass):

    @staticmethod
    def __classcall_private__(cls, *args, **opts):
        r"""
        Ensure that triangular Dyck paths created by the enumerated sets and
        directly are the same and that they are instances of
        :class:`TriangularDyckPath`.

        TESTS::

            sage: tdp = TriangularDyckPath([3,1],[2])
            sage: tdp.parent()
            Triangular Dyck Paths
            sage: tdp2 = TriangularDyckPaths()([3,1],[2])
            sage: tdp2.parent() is tdp.parent()
            True
            sage: type(tdp) is type(tdp2)
            True
        """
        P = TriangularDyckPaths_all()
        return P.element_class(P, *args, **opts)

    def __init__(self, parent, t, p):
        r"""
        TESTS::

            sage: TriangularDyckPath([3,1],[2]).parent()
            Triangular Dyck Paths
        """
        Element.__init__(self, parent)
        self._partition = TriangularPartition(t)
        self._path = Partition(p)
        self._skewpart = SkewPartition([self._partition.partition(), self._path])

    def _repr_(self):
        r"""
        TESTS::

            sage: TriangularDyckPath([3,1],[2])
            [[3, 1], [2]]
        """
        return str(self._skewpart)

    def partition(self):
        r"""
        EXAMPLES::

            sage: TriangularDyckPath([3,1],[2]).partition()
            [3, 1]
        """
        return self._partition

    def path(self):
        r"""
        EXAMPLES::

            sage: TriangularDyckPath([3,1],[2]).path()
            [2]
        """
        return self._path

    def corners(self):
        r"""
        EXAMPLES::

            sage: TriangularDyckPath([3,1],[2]).corners()
            [(0, 1)]
        """
        return self._path.corners()

    def skew_partition(self):
        r"""
        EXAMPLES::

            sage: TriangularDyckPath([3,1],[2]).skew_partition()
            [3, 1] / [2]
        """
        return self._skewpart

    def conjugate(self):
        r"""
        EXAMPLES::

            sage: TriangularDyckPath([3,1],[2]).conjugate()
            [[2, 1, 1], [1, 1]]
        """
        return TriangularDyckPath(self.partition().conjugate(),self._path.conjugate())

    def contains(self, other):
        r"""
        EXAMPLES::

            sage: tdp = TriangularDyckPath([3,1],[2])
            sage: tdp2 = TriangularDyckPath([3,1],[1])
            sage: tdp.contains(tdp2)
            True
        """
        return self.path().contains(other.path())

    def _latex_(self):
        r"""
        TESTS::

            sage: s = latex(TriangularDyckPath([3,1],[2]))
        """
        #return latex(self._skewpart)
        latex.add_package_to_preamble_if_available("tikz")
        return self._grid_latex_()

    def _grid_latex_(self, with_deficit = False, tableau = None):
        partition = self.partition()
        path = self.path()
        maxx = max(partition)
        maxy = partition.length()
        y = 0

        path1 = f"(0,{y})"
        for p in partition:
            path1+= f"--({p},{y})"
            y+=1
            path1+= f"--({p},{y})"
        path1+=f"--(0,{y})--cycle;\n"

        y = 0
        path2 = f"({maxx},0)"
        path3 = "(0,0)"
        for p in path:
            path2+= f"--({p},{y})"
            path3+= f"--({p},{y})"
            y+=1
            path2+= f"--({p},{y})"
            path3+= f"--({p},{y})"
        path2+=f"--(0,{y})--(0,{maxy});\n"
        path3+=f"--(0,{y})--cycle;\n"

        color = "green!10" if with_deficit else "gray!10"

        s = "\\begin{tikzpicture}\n"
        s+= f"\\draw[ultra thick, fill = red!30] " + path1
        s+= f"\\draw[fill = {color}] " + path3
        if with_deficit:
            for y,x in self.deficit_cells(tableau):
                s+=f"\\draw[fill = yellow!50] ({x},{y}) rectangle ({x+1},{y+1});\n"
                s+=f"\\draw[pattern=crosshatch dots, pattern color=gray!50] ({x},{y}) rectangle ({x+1},{y+1});\n"
        if tableau != None:
            for y,x in tableau.cells():
                s+= f"\\node at ({x + .5}, {y +.5}) (int) "
                s+= "{"
                s+= f"${tableau[y][x]}$"
                s+= "};\n"
        s+= f"\\draw[line width = 5px, red] " + path2
        s+= f"\\draw[step=1cm,gray,very thin] (0,0) grid ({maxx},{maxy});\n"
        s+= "\\end{tikzpicture}"
        return s


    def __eq__(self, other):
        return isinstance(other,TriangularDyckPath) and other._skewpart == self._skewpart

    def __hash__(self):
        return hash(self._skewpart)

    def _richcmp_(self, other, op):
        if not isinstance(other, TriangularDyckPath):
            return NotImplemented
        return richcmp(self._skewpart, other._skewpart,op)

    def area(self):
        r"""
        EXAMPLES::

            sage: TriangularDyckPath([3,1],[2]).area()
            2
        """
        return self._skewpart.size()

    def sim(self, tableau = None):
        r"""
        EXAMPLES::

            sage: TriangularDyckPath([3,1],[2]).sim()
            2
            sage: TriangularDyckPath([3,1],[1,1]).sim()
            1
        """
        return len(list(self.similar_cells(tableau)))

    def deficit(self, tableau = None):
        r"""
        EXAMPLES::

            sage: tp = TriangularPartition([3,1])
            sage: tdp = TriangularDyckPath([3,1],[3])
            sage: tdp.deficit()
            1
            sage: tdp.deficit(tableau = tp.top_down_standard_tableau())
            0
        """
        return len(list(self.deficit_cells(tableau)))

    def compose(self, cell, other):
        y, x = cell
        p = self._partition.compose(cell, other.partition())
        p2 = list(other._path)
        p2.extend([0] * (y - len(p2)))
        path = Partition([v + x + 1 for v in p2] + [x+1] + list(self._path))
        return TriangularDyckPath(p, path)

    def path_tamari_rotate(self, line):
        r"""
        EXAMPLES::

            sage: tdp = TriangularDyckPath([3,1],[2])
            sage: tdp.path_tamari_rotate(0)
            [[3, 1], [1]]
        """
        L = self.skew_partition().row_lengths()
        p = list(self.path())
        v = L[line]
        p[line]-=1
        i = line-1
        while i >= 0 and L[i] > v:
            p[i]-=1
            i-=1
        return TriangularDyckPath(self.partition(), p)


    def path_tamari_up(self):
        r"""
        EXAMPLES::

            sage: list(TriangularDyckPath([3,1],[2]).path_tamari_up())
            [[[3, 1], [1]]]
        """
        for l,c in self.path().corners():
            yield self.path_tamari_rotate(l)

    # def triangular_tamari_rotate(self,i,j):
        # if j > 0 and self.is_horizontal_deficit(i,j-1):
            # return self.conjugate().triangular_tamari_rotate(j,i).conjugate()
        # p = list(self.path())
        # p[i]-=1
        # i-=1
        # while i >= 0 and self.is_vertical_deficit(i,j):
            # p[i]-=1
            # i-=1
        # return TriangularDyckPath(self.partition(),p)

    def immediate_deficit_tamari_rotate(self,i,j):
        if j > 0 and self.is_horizontal_deficit(i,j-1):
            return self.conjugate().immediate_deficit_tamari_rotate(j,i).conjugate()
        p = list(self.path())
        L = self.partition().triangular_standard_labels()
        orig = L[(i,j)]
        v = L[(i,j)]
        p[i]-=1
        i-=1
        while i >= 0 and p[i] < self.partition()[i] and v > L[(i,p[i])]:
            v = max(L[(i,j)] for j in range(p[i],self.partition()[i]) if L[(i,j)] < orig)
            p[i]-=1
            i-=1
        return TriangularDyckPath(self.partition(),p)

    def triangular_tamari_rotate(self,i,j):
        if any(self.is_significant_horizontal_deficit(i,j2) for j2 in range(j)):
            return self.conjugate().triangular_tamari_rotate(j,i).conjugate()

        L = self.partition().triangular_standard_labels()
        pinit = list(self.path())
        p = list(self.path())
        p[i]-=1
        v = L[(i,j)]
        p = Partition(p)
        imax = i
        for i2 in range(i-1,-1,-1):
            if self.is_significant_vertical_deficit(i2,j):
                if v > L[(i2,pinit[i2])]:
                    v = max(L[(i2,j2)] for j2 in range(pinit[i2],self.partition()[i2]) if L[(i2,j2)] < v)
                    c = max((c for c in p.corners() if c[0] >= i2 and c[0] < imax and c[1] >= j), key = lambda c: L[c])
                    p = list(p)
                    p[c[0]] -= 1
                    p = Partition(p)
                    imax = c[0]
                else:
                    break

        return TriangularDyckPath(self.partition(),p)

    def slope_tamari_rotate(self,i,j):
        p = list(self.path())
        v = self.cell_mean_slope(i,j)
        p[i]-=1
        i = i-1
        while i >= 0 and self.cell_mean_slope(i, j) > v:
            #v = self.cell_mean_slope(i, j)
            p[i]-=1
            i-=1
        return TriangularDyckPath(self.partition(), p)

    @cached_method
    def triangular_tamari_removable_cells(self):
        if self.path().size() > 0:
            L = self.partition().triangular_standard_labels()
            cell = max(self.corners(), key = lambda c:L[c])
            s = self.triangular_tamari_rotate(*cell).triangular_tamari_removable_cells()
            s.add(cell)
            return s
        return set()

    @cached_method
    def triangular_tamari_plus_cells(self):
        rcells = self.triangular_tamari_removable_cells()
        return {c for c in self.similar_cells() if c[0]==0 and not c in rcells}

    @cached_method
    def triangular_tamari_minus_cells(self):
        s = set()
        rcells = self.triangular_tamari_removable_cells()
        L = self.partition().triangular_standard_labels()
        path = self.path()
        conj = self.path().conjugate()
        for i,j in self.skew_partition().cells():
            for i2 in range(i):
                if (i2,j) in rcells and L[(i,j)] < L[(i2,path[i2]-1)]:
                    s.add((i,j))
                    break
            else:
                for j2 in range(j):
                    if (i,j2) in rcells and L[(i,j)] < L[(conj[j2]-1,j2)]:
                        s.add((i,j))
                        break
        return s



    # ~ def triangular_tamari_rotate(self,i,j):
        # ~ if j > 0 and self.is_horizontal_deficit(i,j-1):
            # ~ return self.conjugate().triangular_tamari_rotate(j,i).conjugate()
        # ~ p = list(self.path())
        # ~ L = self.partition().triangular_standard_labels()
        # ~ v = L[(i,j)]
        # ~ p[i]-=1
        # ~ i-=1
        # ~ while i >= 0 and p[i] < self.partition()[i] and v > L[(i,p[i])]:
            # ~ v = max(L[(i,j)] for j in range(p[i],self.partition()[i]) if L[(i,j)] < v)
            # ~ p[i]-=1
            # ~ i-=1
        # ~ return TriangularDyckPath(self.partition(),p)

    # def triangular_tamari_rotate(self,i,j):
        # if j > 0 and self.is_horizontal_deficit(i,j-1):
            # return self.conjugate().triangular_tamari_rotate(j,i).conjugate()
        # p = list(self.path())
        # L = self.partition().triangular_standard_labels()
        # orig = L[(i,j)]
        # v = orig
        # p[i]-=1
        # i-=1
        # while i >= 0 and p[i] < self.partition()[i]:
            # w = max((L[(i,j)] for j in range(p[i],self.partition()[i]) if L[(i,j)] < orig), default = orig+1)
            # if w > v:
                # break
            # v = w
            # p[i]-=1
            # i-=1
        # return TriangularDyckPath(self.partition(),p)

    def triangular_tamari_up(self):
        for i,j in self.corners():
            yield self.triangular_tamari_rotate(i,j)

    def immediate_deficit_tamari_up(self):
        for l,c in self.corners():
            yield self.immediate_deficit_tamari_rotate(l,c)

    def slope_tamari_up(self):
        for l,c in self.corners():
            yield self.slope_tamari_rotate(l,c)

    # def triangular_schur_rotate(self,i,j):
        # L = self.partition().triangular_standard_labels()
        # v = L[(i,j)]
        # d = set()
        # d.update((i,j2) for j2 in range(j) if self.is_horizontal_deficit(i,j2))
        # d.update((i2,j) for i2 in range(i) if self.is_vertical_deficit(i2,j))
        # prev = len(d)
        # nd = len(d)
        # diff = 0
        # p = list(self.path())
        # p[i]-=1
        # tdp = TriangularDyckPath(self.partition(),p)
        # while nd > 0:
            # for k in range(nd):
                # corners = tdp.corners()
                # corners.sort(key = lambda c: L[c])
                # corner = corners[-(diff+1)]
                # d.update((corner[0],j2) for j2 in range(corner[1]) if self.is_horizontal_deficit(corner[0],j2))
                # d.update((i2,corner[1]) for i2 in range(corner[0]) if self.is_vertical_deficit(i2,corner[1]))
                # if L[corner] > v:
                    # diff+=1
                # p = list(tdp.path())
                # p[corner[0]]-=1
                # tdp = TriangularDyckPath(self.partition(),p)
            # nd = len(d) - prev
            # prev = len(d)
        # return tdp

    def triangular_schur_rotate(self,i,j):
        L = self.partition().triangular_standard_labels()
        v = L[(i,j)]
        d = set()
        d.update((i,j2) for j2 in range(j) if self.is_horizontal_deficit(i,j2))
        d.update((i2,j) for i2 in range(i) if self.is_vertical_deficit(i2,j))
        prev = len(d)
        nd = len(d)
        p = list(self.path())
        p[i]-=1
        tdp = TriangularDyckPath(self.partition(),p)
        while nd > 0:
            for k in range(nd):
                corners = [c for c in tdp.corners() if L[c] < v]
                corners.sort(key = lambda c: L[c])
                if len(corners) == 0:
                    return tdp
                corner = corners[-1]
                d.update((corner[0],j2) for j2 in range(corner[1]) if self.is_horizontal_deficit(corner[0],j2))
                d.update((i2,corner[1]) for i2 in range(corner[0]) if self.is_vertical_deficit(i2,corner[1]))
                p = list(tdp.path())
                p[corner[0]]-=1
                tdp = TriangularDyckPath(self.partition(),p)
            nd = len(d) - prev
            prev = len(d)
        return tdp

    def triangular_schur_up(self):
        for i,j in self.corners():
            yield self.triangular_schur_rotate(i,j)

    def cell_min_slope(self,i,j):
        r"""
        EXAMPLES::

            sage: tdp = TriangularDyckPath([3,1],[2])
            sage: tdp.cell_min_slope(0,0)
            0
        """
        p = self._path
        return p.leg_length(i,j) / p.hook_length(i,j)

    def cell_max_slope(self,i,j):
        r"""
        EXAMPLES::

            sage: tdp = TriangularDyckPath([3,1],[2])
            sage: tdp.cell_max_slope(0,0)
            1/2
        """
        p = self._path
        return (p.leg_length(i,j)+1) / p.hook_length(i,j)

    def cell_mean_slope(self,i,j):
        r"""
        EXAMPLES::

            sage: tdp = TriangularDyckPath([3,1],[2])
            sage: tdp.cell_mean_slope(0,0)
            1/4
        """
        return (self.cell_min_slope(i,j) + self.cell_max_slope(i,j))/2

    def is_sim_cell(self,i,j, tableau = None):
        r"""
        EXAMPLES::

            sage: tdp = TriangularDyckPath([3,1],[2])
            sage: tdp.is_sim_cell(0,0)
            True
            sage: tdp = TriangularDyckPath([3,1],[1,1])
            sage: tdp.is_sim_cell(0,0)
            False
        """
        if tableau is None:
            tau = self.partition()
            return self.cell_min_slope(i,j) < tau.mean_slope() and tau.mean_slope() <= self.cell_max_slope(i,j)
        else:
            return not self.is_deficit_cell(i,j,tableau)

    def is_deficit_cell(self,i,j, tableau = None):
        r"""
        EXAMPLES::

            sage: tdp = TriangularDyckPath([3,1],[3])
            sage: tdp.is_deficit_cell(0,0)
            True
            sage: tdp.is_deficit_cell(0,0, tableau = tdp.partition().top_down_standard_tableau())
            False
        """
        col_ext = [c for c in self.skew_partition().cells() if c[1] == j]
        row_ext = [c for c in self.skew_partition().cells() if c[0] == i]
        col_int = [c for c in self.path().cells() if c[1] == j]
        row_int = [c for c in self.path().cells() if c[0] == i]
        if tableau is None:
            L = self.partition().triangular_standard_tableau()
        else:
            L = tableau
        return any(L.entry(c1) < L.entry(c2) for c1 in col_ext for c2 in row_int) or any(L.entry(c1) < L.entry(c2) for c1 in row_ext for c2 in col_int)

    # def is_vertical_lead_deficit(self,i,j):
        # return self.is_vertical_deficit(i,j) and all(not self.is_vertical_deficit(i,j2) for j2 in range(j+1,self.path()[i]))

    # def is_horizontal_lead_deficit(self,i,j):
        # return self.is_horizontal_deficit(i,j) and all(not self.is_horizontal_deficit(i2,j) for i2 in range(i+1,self.path().conjugate()[j]))

    # def is_lead_deficit_cell(self,i,j):
        # return self.is_vertical_lead_deficit(i,j) or self.is_horizontal_lead_deficit(i,j)

    def is_horizontal_deficit(self,i,j):
        leg = [c for c in self.skew_partition().cells() if c[1] == j]
        arm = self.path().arm_cells(i,j)
        if len(arm) == 0:
            return False
        corner = arm[-1]
        L = self.partition().triangular_standard_labels()
        return any(L[c] < L[corner] for c in leg)

    def is_vertical_deficit(self,i,j):
        arm = [c for c in self.skew_partition().cells() if c[0] == i]
        leg = self.path().leg_cells(i,j)
        if len(leg) == 0:
            return False
        corner = leg[-1]
        L = self.partition().triangular_standard_labels()
        return any(L[c] < L[corner] for c in arm)

    def is_significant_horizontal_deficit(self,i,j):
        leg = [c for c in self.skew_partition().cells() if c[1] == j]
        arm = self.path().arm_cells(i,j)
        if len(arm) == 0:
            return False
        corner = arm[-1]
        L = self.partition().triangular_standard_labels()
        for c in leg:
            if L[c] < L[corner]:
                for c2 in self.partition().horizontal_rectangle_cells(c,corner):
                    if c2[0] == corner[0] + 1 and L[c2] > L[corner]:
                        break
                    if L[c2] > L[c] and not self.is_label_in_path(L[c2]):
                        break
                else:
                    return True
        return False

    def is_significant_vertical_deficit(self, i, j):
        arm = [c for c in self.skew_partition().cells() if c[0] == i]
        leg = self.path().leg_cells(i,j)
        if len(leg) == 0:
            return False
        corner = leg[-1]
        L = self.partition().triangular_standard_labels()
        for c in arm:
            if L[c] < L[corner]:
                for c2 in self.partition().vertical_rectangle_cells(corner,c):
                    if c2[1] == corner[1] + 1 and L[c2] > L[corner]:
                        break
                    if L[c2] > L[c] and not self.is_label_in_path(L[c2]):
                        break
                else:
                    return True
        return False

    def path_labels(self):
        L = self.partition().triangular_standard_labels()
        return set(L[c] for c in self.skew_partition().cells())

    def is_label_in_path(self,v):
        return v in self.path_labels()

    def similar_cells(self, tableau = None):
        r"""
        EXAMPLES::

            sage: tdp = TriangularDyckPath([3,1],[3])
            sage: list(tdp.similar_cells())
            [(0, 1), (0, 2)]
        """
        for i,j in self._path.cells():
            if self.is_sim_cell(i,j, tableau):
                yield (i,j)

    def deficit_cells(self, tableau = None):
        r"""
        EXAMPLES::

            sage: tdp = TriangularDyckPath([3,1],[3])
            sage: list(tdp.deficit_cells())
            [(0, 0)]
            sage: list(tdp.deficit_cells(tableau = tdp.partition().top_down_standard_tableau()))
            []
        """
        for i,j in self._path.cells():
            if self.is_deficit_cell(i,j, tableau):
                yield (i,j)

    # def lead_deficit_cells(self):
        # for i,j in self._path.cells():
            # if self.is_lead_deficit_cell(i,j):
                # yield (i,j)

    def path_lattice_interval(self, other):
        r"""
        EXAMPLES::

            sage: tdp1 = TriangularDyckPath([3,1],[3])
            sage: tdp2 = TriangularDyckPath([3,1],[1])
            sage: tdp1.path_lattice_interval(tdp2)
            Finite lattice containing 3 elements
        """
        L = self.partition().path_tamari_lattice()
        if not L.le(self,other):
            return
        return LatticePoset(L.subposet(L.interval(self,other)))


    def path_distance(self,other):
        r"""
        EXAMPLES::

            sage: tdp1 = TriangularDyckPath([3,1],[3])
            sage: tdp2 = TriangularDyckPath([3,1],[1])
            sage: tdp1.path_distance(tdp2)
            2
        """
        I = self.path_lattice_interval(other)
        if I is None:
            return -1
        return max(len(c) for c in I.maximal_chains()) - 1

    def triangular_poset_interval(self, other):
        L = self.partition().triangular_tamari_poset()
        if not L.le(self,other):
            return
        return L.subposet(L.interval(self,other))

    def triangular_distance(self,other):
        I = self.triangular_poset_interval(other)
        if I is None:
            return -1
        return max(len(c) for c in I.maximal_chains()) - 1

    def descents_composition(self):
        r"""
        EXAMPLES::

            sage: tdp = TriangularDyckPath([3,1],[3])
            sage: tdp.descents_composition()
            [2, 1]
            sage: TriangularDyckPath([3,1],[1,1]).descents_composition()
            [1, 2]
        """
        p = list(self.path()) + [0] * (self.partition().length() - self.path().length() + 1)
        c = []
        p.reverse()
        v = 1
        for i in range(1,len(p)):
            if p[i] == p[i-1]:
                v+=1
            else:
                c.append(v)
                v = 1
        c.append(v)
        return Composition(c)

    def descents_partition(self):
        r"""
        EXAMPLES::

            sage: tdp = TriangularDyckPath([3,1],[3])
            sage: tdp.descents_partition()
            [2, 1]
            sage: TriangularDyckPath([3,1],[1,1]).descents_partition()
            [2, 1]
        """
        return Partition(reversed(sorted(self.descents_composition())))

    def cells_RSK(self):
        L = self.partition().triangular_standard_labels()
        word1 = self.skew_partition().cells()
        word1.sort(key = lambda x: -L[x])
        word2 = [c for c in self.path().cells() if not self.is_deficit_cell(*c)]
        word2.sort(key = lambda x : -L[x])
        word = word1 + word2
        return RSK([self.partition().size() + 1 - L[x] for x in word])

    def cells_RSK_partition(self):
        return self.cells_RSK()[0].shape()



class TriangularDyckPaths(UniqueRepresentation, Parent):

    def __call__(self, *args, **keywords):

        if isinstance(args[0], TriangularDyckPath):
            return args[0]

        return super(TriangularDyckPaths, self).__call__(*args, **keywords)

    @staticmethod
    def __classcall_private__(cls, t=None):
        r"""
        TESTS::

            sage: TriangularDyckPaths()
            Triangular Dyck Paths
            sage: TriangularDyckPaths([3,1])
            Triangular Dyck Paths of [3, 1]
        """
        if t is None:
            return TriangularDyckPaths_all()

        return TriangularDyckPaths_partition(TriangularPartition(t))

class TriangularDyckPaths_all(DisjointUnionEnumeratedSets, TriangularDyckPaths):

    def __init__(self):
        DisjointUnionEnumeratedSets.__init__(
            self, Family(TriangularPartitions(), TriangularDyckPaths_partition),
            facade=True, keepkey=False,
            category=EnumeratedSets())

    def _repr_(self):
        r"""
        TESTS::

            sage: TriangularDyckPaths()
            Triangular Dyck Paths
        """
        return "Triangular Dyck Paths"

    def _element_constructor_(self, t, p):
        r"""
        EXAMPLES::

            sage: TDP = TriangularDyckPaths()
            sage: TDP([3,1],[2])
            [[3, 1], [2]]
        """
        return self.element_class(self, t, p)

    def __contains__(self, x):
        r"""
        TESTS::

            sage: TDP = TriangularDyckPaths()
            sage: 1 in TDP
            False
            sage: TDP([3,1],[2]) in TDP
            True
        """
        return isinstance(x, self.element_class)

    Element = TriangularDyckPath

class TriangularDyckPaths_partition(TriangularDyckPaths):
    r"""
    The enumerated set of triangular Dyck paths for a given partition
    """
    def __init__(self, partition):
        # sage: for tp in TriangularPartitions(6).some_elements(): TestSuite(TriangularDyckPaths(tp)).run()
        # fail because of pickling
        r"""
        TESTS::

            sage: TDP = TriangularDyckPaths([3,1])
            sage: assert TDP is TriangularDyckPaths([3,1])
        """
        super(TriangularDyckPaths_partition, self).__init__(category=FiniteEnumeratedSets())

        self._partition = partition

    def partition(self):
        r"""
        EXAMPLES::

            sage: TDP = TriangularDyckPaths([3,1])
            sage: TDP.partition()
            [3, 1]
        """
        return self._partition

    def _repr_(self):
        r"""
        TESTS::

            sage: TriangularDyckPaths([3,1])
            Triangular Dyck Paths of [3, 1]
        """
        return "Triangular Dyck Paths of {}".format(self._partition)

    def __contains__(self, x):
        r"""
        TESTS::

            sage: TDP = TriangularDyckPaths([3,1])
            sage: TriangularDyckPath([3,1],[2]) in TDP
            True
            sage: TriangularDyckPath([4,1],[2]) in TDP
            False
            sage: 1 in TDP
            False
        """
        return isinstance(x, self.element_class) and x.partition() == self._partition

    @lazy_attribute
    def _parent_for(self):
        r"""
        The parent of the element generated by ``self``.

        TESTS::

            sage: TDP = TriangularDyckPaths([3,1])
            sage: TDP._parent_for
            Triangular Dyck Paths
        """
        return TriangularDyckPaths_all()

    # This is needed because this is a facade parent via DisjointUnionEnumeratedSets
    @lazy_attribute
    def element_class(self):
        r"""
        TESTS::

            sage: TDP = TriangularDyckPaths([3,1])
            sage: TDP.element_class
            <class '__main__.TriangularDyckPaths_all_with_category.element_class'>
            sage: TDP.first().__class__ == TriangularDyckPaths().first().__class__
            True
        """
        return self._parent_for.element_class

    def _element_constructor_(self, path):
        r"""
        EXAMPLES::

            sage: TDP = TriangularDyckPaths([3,1])
            sage: TDP([2])
            [[3, 1], [2]]
        """
        return self.element_class(self, self._partition, path)

    def __iter__(self):
        r"""
        EXAMPLES::

            sage: list(TriangularDyckPaths([2]))
            [[[2], []], [[2], [1]], [[2], [2]]]
        """
        if self.partition().size() == 0:
            yield TriangularDyckPath(self.partition(), [])
        else:
            for p in TriangularDyckPaths(self.partition().interior()):
                yield TriangularDyckPath(self.partition(), p.path())

            for skewalpha, beta, cell in self.partition().diagonal_decompositions():
                for p1 in TriangularDyckPaths(skewalpha[1]):
                    for p2 in TriangularDyckPaths(beta):
                        yield TriangularDyckPath(skewalpha[0], p1.path()).compose(cell, p2)


