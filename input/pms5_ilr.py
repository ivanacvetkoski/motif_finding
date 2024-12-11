import numpy as np
from pulp import LpMaximize, LpProblem, LpVariable, lpSum, LpStatus


def solve_ilp(n1, n2, n3, n4, n5, p1, p2, p3):
    prob = LpProblem("ILP", LpMaximize)

    N1a = LpVariable('N1a', lowBound=0, cat='Integer')
    N2a = LpVariable('N2a', lowBound=0, cat='Integer')
    N2b = LpVariable('N2b', lowBound=0, cat='Integer')
    N3a = LpVariable('N3a', lowBound=0, cat='Integer')
    N3b = LpVariable('N3b', lowBound=0, cat='Integer')
    N4a = LpVariable('N4a', lowBound=0, cat='Integer')
    N4b = LpVariable('N4b', lowBound=0, cat='Integer')
    N5a = LpVariable('N5a', lowBound=0, cat='Integer')
    N5b = LpVariable('N5b', lowBound=0, cat='Integer')
    N5c = LpVariable('N5c', lowBound=0, cat='Integer')

    prob += n1-N1a + n2-N2a + n3-N3a + n4-N4b + n5-N5a <= p1
    prob += n1-N1a + n2-N2a + n3-N3b + n4-N4a + n5-N5b <= p2
    prob += n1-N1a + n2-N2b + n3-N3a + n4-N4a + n5-N5c <= p3
    prob += N1a <= n1
    prob += N2a + N2b <= n2
    prob += N3a + N3b <= n3
    prob += N4a + N4b <= n4
    prob += N5a + N5b + N5c <= n5

    prob += lpSum([N1a, N2a, N2b, N3a, N3b, N4a, N4b, N5a, N5b, N5c])

    prob.solve()

    return LpStatus[prob.status] == 'Optimal'


def generate_8d_table(l, d):
    shape = (l + 1, l + 1, l + 1, l + 1, l + 1, d + 1, d + 1, d + 1)
    table = np.zeros(shape, dtype=bool)

    with open(f"ilr_{l}_{d}.txt", "w") as file:
        for n1 in range(l + 1):
            for n2 in range(l + 1):
                for n3 in range(l + 1):
                    for n4 in range(l + 1):
                        for n5 in range(l + 1):
                            if n1 + n2 + n3 + n4 + n5 > l:
                                continue
                            for p1 in range(d + 1):
                                for p2 in range(d + 1):
                                    for p3 in range(d + 1):
                                        if solve_ilp(n1, n2, n3, n4, n5, p1, p2, p3):
                                            file.write(f"{n1} {n2} {n3} {n4} {n5} {p1} {p2} {p3} true\n")
                                        else:
                                            file.write(f"{n1} {n2} {n3} {n4} {n5} {p1} {p2} {p3} false\n")

    return table

generate_8d_table(8, 1)