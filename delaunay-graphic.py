import tkinter as tk
from itertools import chain
from math import sqrt

from delaunay import *

fenetre = tk.Tk()
c = tk.Canvas(fenetre, height=800, width=800, bg="white")
c.pack()

points = []

def ajouter_point(event):
    points.append((event.x, event.y))
    c.create_oval(event.x-3, event.y-3, event.x+3, event.y+3, fill="black")

c.bind("<Button-1>", ajouter_point)

def reset():
    global points
    c.delete(tk.ALL)
    points = []

def tracer(event=None):
    """Trace la triangulation de Delaunay de l'ensemble des points."""
    succ, _ = delaunay(points)
    for (a, b) in succ:
        c.create_line(*a, *b, width=1.5)

def tracer_cercles(points):
    for x, y in points:
        c.create_oval(x-3, y-3, x+3, y+3, fill="black")
    aretes = delaunay(points)
    for (a, b) in aretes:
        c.create_line(*a, *b, width=1.5)
        d = aretes[a, b]
        if (b, d) in aretes:
            (xo, yo), r = cercle_circonscrit(a, b, d)
            c.create_oval(xo - r, yo - r, xo + r, yo + r, width=1)    

def cercle_circonscrit(a, b, c):
    """Renvoie le couple (o, r) avec o le centre du cercle circonscrit à
    a, b, c et r le rayon."""
    xa, ya = a
    xb, yb = b
    xc, yc = c
    d = 2 * (xa * (yb - yc) + xb * (yc - ya) + xc * (ya - yb))
    xo = ((xa**2 + ya**2)*(yb - yc) + (xb**2 + yb**2)*(yc - ya) + (xc**2 + yc**2)*(ya - yb))/d   
    yo = ((xa**2 + ya**2)*(xc - xb) + (xb**2 + yb**2)*(xa - xc) + (xc**2 + yc**2)*(xb - xa))/d
    r = sqrt((xa - xo) ** 2 + (ya - yo) ** 2)
    return (xo, yo), r

def test_delaunay(points):
    succ = delaunay(points)
    for (a, b), c in succ.items():
        if (b, c) in succ:
            for d in points:
                if d not in (a, b, c) and orientation(a, b, c) == DIRECT:
                    assert position_cercle_circonscrit(a, b, c, d) != DEDANS, \
                           (a, b, c, d)

def step_by_step_delaunay(points):
    """Calcule la triangulation de Delaunay d'une liste de points distincts du
    plan.

    Utilise l'algorithme 'divide and conquer' de Lee et Schachter :
    http://www.personal.psu.edu/cxc11/AERSP560/DELAUNEY/13_Two_algorithms_Delauney.pdf"""
    succ = {}
    pred = {}
    first = {}

    def delete(a, b):
        sa = succ.pop((a, b))
        sb = succ.pop((b, a))
        pa = pred.pop((a, b))
        pb = pred.pop((b, a))
        succ[a, pa] = sa
        succ[b, pb] = sb
        pred[a, sa] = pa
        pred[b, sb] = pb
        return ("delete", (a, b))

    def insere(a, b, sa, pb):
        """Insère l'arête (a, b), avec sa et pb les points tels qu'à la fin de
        l'opération, succ[a, b] = sa et pred[b, a] = pb."""
        pa = pred[a, sa]
        sb = succ[b, pb]
        succ[a, pa] = b
        succ[a, b] = sa
        pred[a, sa] = b
        pred[a, b] = pa
        pred[b, sb] = a
        pred[b, a] = pb
        succ[b, pb] = a
        succ[b, a] = sb

    def insere_line(a, b):
        return ("insert", (a, b))

    def common_tangent(x0, y0):
        """Si nous avons calculé la tangente commune à deux ensembles de points
        X et Y séparés par une droite, et que x0 et y0 sont les points de X et Y
        les plus proches de la droite, calcule la tangente commune des
        enveloppes convexes de X et Y passant par les points x de X et y de Y
        et telle que tous les points de X et de Y se trouvent à gauche du
        segment orienté (x, y)."""
        x, y = x0, y0
        z0 = first[y]
        z1 = first[x]
        assert (x, z1) in succ, (x, z1, succ)
        assert (z1, x) in succ, (z1, x, succ)
        z2 = pred[x, z1]
        while True:
            if orientation(x, y, z0) == INDIRECT:
                y, z0 = z0, succ[z0, y]
            elif orientation(x, y, z2) == INDIRECT:
                x, z2 = z2, pred[z2, x]
            else:
                return (x, y)

    def merge(x, y):
        """En prenant les mêmes notations que la docstring précédente, on
        fusionne les triangulations de Delaunay des ensembles X et Y.
        Pour cela, on met à jour les informations à propos de succ, pred et
        first."""
        insere(x, y, first[x], pred[y, first[y]])
        yield "insert", (x, y)
        first[x] = y
        
        while True: 
            if orientation(x, y, pred[y, x]) == DIRECT:
                y1 = pred[y, x]
                y2 = pred[y, y1]
                while position_cercle_circonscrit(x, y, y1, y2) == DEDANS:
                    yield delete(y, y1)
                    y1 = y2
                    y2 = pred[y, y1]
            else:
                y1 = None


            if orientation(x, y, succ[x, y]) == DIRECT:
                x1 = succ[x, y]
                x2 = succ[x, x1]
                while position_cercle_circonscrit(x, y, x1, x2) == DEDANS:
                    yield delete(x, x1)
                    x1 = x2
                    x2 = succ[x, x1]
            else:
                x1 = None

            if x1 is None and y1 is None:
                break
            elif x1 is None:
                #on trace l'arête (x, y1)
                insere(y1, x, y, y)
                yield "cercle", (x, y1, y)
                y = y1
            elif y1 is None:
                #on trace l'arête (y, x1)
                insere(y, x1, x, x)
                yield "cercle", (y, x1, x)
                x = x1
            elif position_cercle_circonscrit(x, y, y1, x1) == DEDANS:
                insere(y, x1, x, x)
                yield "cercle", (y, x1, x)
                x = x1
            else:
                insere(y1, x, y, y)
                yield "cercle", (x, y1, y)
                y = y1

        #Finalement, on met à jour les informations sur l'enveloppe convexe
        #commune.
        first[y] = x

    def compute(points):
        """Calcule la triangulation de Delaunay et l'enveloppe convexe d'une
        liste de points en mettant à jour les structures de données
        succ, pred et first.

        Préconditions :
        -la liste points contient au moins deux éléments
        -les points sont tous distincts
        -les points sont orientés dans l'ordre lexicographique"""
        n = len(points)
        
        if n == 2:
            [a, b] = points
            #S'il n'y a que deux points, on trace le segment [a, b].
            succ[a, b] = pred[a, b] = b
            succ[b, a] = pred[b, a] = a
            first[a] = b
            first[b] = a
            yield "insert", (a, b)
            
        elif n == 3:
            [a, b, c] = points
            #Trois cas : le triangle (a, b, c) est direct
            #            le triangle (a, b, c) est indirect
            #            les points (a, b, c) sont alignés
            if orientation(a, b, c) == DIRECT:
                succ[a, c] = succ[c, a] = pred[a, c] = pred[c, a] = b
                succ[a, b] = succ[b, a] = pred[a, b] = pred[b, a] = c
                succ[b, c] = succ[c, b] = pred[b, c] = pred[c, b] = a
                first[a] = b
                first[b] = c
                first[c] = a
                yield "insert", (a, b)
                yield "insert", (a, c)
                yield "insert", (b, c)
            elif orientation(a, b, c) == INDIRECT:
                succ[a, b] = succ[b, a] = pred[a, b] = pred[b, a] = c
                succ[a, c] = succ[c, a] = pred[a, c] = pred[c, a] = b
                succ[b, c] = succ[c, b] = pred[b, c] = pred[c, b] = a
                first[a] = c
                first[b] = a
                first[c] = b
                yield "insert", (a, b)
                yield "insert", (a, c)
                yield "insert", (b, c)
            else:
                #on trie de façon à ce que b soit le point du milieu et a et
                #c les extrêmes.
                [a, b, c] = sorted(points)
                #Si a, b, c sont alignés, on ne trace que les arêtes (a, b)
                #et (a, c). Quant à l'enveloppe convexe, on ne considère que
                #les points extrêmes a et c.
                succ[a, b] = pred[a, b] = succ[c, b] = pred[c, b] = b
                succ[b, a] = pred[b, a] = c
                succ[b, c] = pred[b, c] = a
                #On considère que first[a] = b et non c puisqu'il n'y a pas
                #d'arête (a, c). Idem pour first[c]. 
                first[a] = b
                first[c] = b
                yield "insert", (a, b)
                yield "insert", (b, c)
                
        else: #S'il y a au moins 4 points :
            var_x, var_y = variance_xy(points)
            if var_y < var_x:
                #On sépare les points selon une droite verticale :
                med = pseudo_mediane(points)
                left = [p for p in points if p < med]
                right = [p for p in points if p >= med]
                for event in compute(left):
                    yield event
                for event in compute(right):
                    yield event
                x, y = common_tangent(max(left), min(right))
                for event in merge(x, y):
                    yield event
            else:
                #On sépare les points selon une droite horizontale :
                med = pseudo_mediane(points, key=inv)
                down = [p for p in points if inv(p) < inv(med)]
                up = [p for p in points if inv(p) >= inv(med)]
                for event in compute(down):
                    yield event
                for event in compute(up):
                    yield event
                x, y = common_tangent(max(down, key=inv), min(up, key=inv))
                for event in merge(x, y):
                    yield event
            
    for event in compute(points):
        yield event


edges = {}
num_cercle = None
points = []

iterateur = None
def iter_del(event=None):
    if iterateur is None:
        print("Placez des points à l'aide de la souris.")
        return
    global num_cercle
    if num_cercle is not None:
        c.delete(num_cercle)
        num_cercle = None
    try:
        act, pts = next(iterateur)
        if act == "insert":
            print(act, pts)
            a, b = pts
            num = c.create_line(a[0], a[1], b[0], b[1], width=1.5)
            edges[a, b] = num
        elif act == "cercle":
            a, b, d = pts
            num = c.create_line(a[0], a[1], b[0], b[1], width=1.5)
            edges[a, b] = num
            (xo, yo), r = cercle_circonscrit(a, b, d)
            num_cercle = c.create_oval(xo-r, yo-r, xo+r, yo+r,
                                       width=1.5)
        else:
            print(act, pts)
            a, b = pts
            if (a, b) in edges:
                num = edges.pop((a, b))
            else:
                num = edges.pop((b, a))
            c.delete(num)
    except StopIteration:
        pass

c.focus_set()
c.bind('<Return>', tracer)

def sbs_button_start(event=None):
    global iterateur
    iterateur = step_by_step_delaunay(points)
    print("start")

b = tk.Button(fenetre, text="reset", command=reset)
b2 = tk.Button(fenetre, text="step by step",
               command=sbs_button_start)
c.bind("<space>", iter_del)

b.pack()
b2.pack()

print("Placez des points à l'aide de la souris.\n"
      "Appuyez sur Entrée pour tracer la triangulation de Delaunay,\n"
      "Ou sur step by step pour voir le déroulement de l'algorithme\n"
      "en appuyant sur la barre d'espace pour chaque étape.")

fenetre.mainloop()
