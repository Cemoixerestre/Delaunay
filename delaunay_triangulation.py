###############################################################################
########################## TRIANGULATION DE DELAUNAY ##########################
###############################################################################

DIRECT = 1
ALIGNES = 0
INDIRECT = -1

DEDANS = 1
CERCLE = 0
DEHORS = -1

def orientation(a, b, c):
    xa, ya = a
    xb, yb = b
    xc, yc = c
    #le déterminant des vecteurs ab et ac :
    d = (xb - xa) * (yc - ya) - (yb - ya) * (xc - xa)
    if d > 0:
        return DIRECT
    elif d == 0:
        return ALIGNES
    return INDIRECT

def position_cercle_circonscrit(a, b, c, d):
    """Renvoie la position du point d par rapport au cercle circonscrit à
    a, b, c. Cette position prend ses valeurs parmi DEDANS, CECLE et DEHORS.

    Précondition : les points a, b, c sont orientés dans le sens direct."""
    ax, ay = a
    bx, by = b
    cx, cy = c
    dx, dy = d
    #voir : https://fr.wikipedia.org/wiki/Triangulation_de_Delaunay#Algorithmes
    m = [[ax - dx, ay - dy, (ax - dx) ** 2 + (ay - dy) ** 2],
         [bx - dx, by - dy, (bx - dx) ** 2 + (by - dy) ** 2],
         [cx - dx, cy - dy, (cx - dx) ** 2 + (cy - dy) ** 2]]
    #le déterminant de la matrice m :
    d = m[0][0]*m[1][1]*m[2][2] + m[0][1]*m[1][2]*m[2][0] + m[0][2]*m[1][0]*m[2][1] -\
        m[0][0]*m[1][2]*m[2][1] - m[0][1]*m[1][0]*m[2][2] - m[0][2]*m[1][1]*m[2][0]
    if d > 0:
        return DEDANS
    elif d == 0:
        return CERCLE
    return DEHORS

def mediane(l, key=None):
    """La médiane d'une liste l."""
    n = len(l)
    l = sorted(l, key=key)
    return l[n//2]

def _pseudo_mediane(l, start, end, key=None, k=7, limmed=100):
    """Calcul d'un médian linéaire d'une liste l.

    Si la taille de l est inférieur à limmed, on calcule la médiane de l.
    Sinon, on partitionne l en k listes, on calcule récursivement les médians
    linéaires des sous-listes et on revoie la 'médiane des médianes'.

    Le calcul se fait en O(n) contre O(n*log(n)) pour la médiane.
    Les valeurs par défaut de k et limmed ont étées déterminées
    empiriquement de façon à optimiser le temps de calcul pour un grand nombre
    de points."""
    if end - start <= limmed:
        return mediane(l[start:end], key)
    sous_medianes = [_pseudo_mediane(l, start+(i*(end-start))//k,
                                        start+((i+1)*(end-start))//k,
                                        key, k, limmed)
                     for i in range(k)]
    return mediane(sous_medianes, key)

def pseudo_mediane(l, key=None, k=7, limmed=100):
    return _pseudo_mediane(l, 0, len(l), key, k, limmed)

def inv(p):
    """Inverse les coordonnées d'un point p.
    Cette fonction est utilisée pour comparer les points par rapport à leurs
    ordonnées et non leurs abscisses, comme avec l'ordre lexicographique sur
    les tuples."""
    x, y = p
    return (y, x)

def variance_xy(points):
    """Calcule la variance des coordonées x et y des points."""
    n = len(points)
    sum_x, sum_y, sum_sqx, sum_sqy = 0, 0, 0, 0
    for (x, y) in points:
        sum_x += x
        sum_y += y
        sum_sqx += x ** 2
        sum_sqy += y ** 2
    var_x = (sum_sqx / n) - (sum_x / n) ** 2
    var_y = (sum_sqy / n) - (sum_y / n) ** 2
    return var_x, var_y

def delaunay_triangulation(points):
    """Calcule la triangulation de Delaunay d'une liste de points distincts du
    plan.

    Utilise l'algorithme 'divide and conquer' de Lee et Schachter :
    http://www.personal.psu.edu/cxc11/AERSP560/DELAUNEY/13_Two_algorithms_Delauney.pdf

    Renvoie deux dictionnaires succ et pred, définis ci-dessous."""

    #succ : un dictionnaire qui associe à des couples de points (a, b) un point
    #c. Plus précisément : si (a, b) est une arête du graphe, si l'on considère
    #la liste des arêtes partant de a, si l'arête (a, c) apparait immédiatement
    #après l'arête (a, b) en tournant dans le sens trigonométrique, alors
    #succ[a, b] = c.
    succ = {}
    #idem que pour succ, mais en tournant dans le sens des aiguilles d'une
    #montr.
    pred = {}
    #Si a est un point de l'enveloppe convexe d'un ensemble de points, si l'on
    #parcourt cette enveloppe dans le sens trigonométrique et que le point b est
    #relié à a dans ce sens, alors first[a] = b.
    first = {}

    def delete(a, b):
        """Supprime l'arête (a, b) et conserve les invariants des structures
        de données pred et succ."""
        sa = succ.pop((a, b))
        sb = succ.pop((b, a))
        pa = pred.pop((a, b))
        pb = pred.pop((b, a))
        succ[a, pa] = sa
        succ[b, pb] = sb
        pred[a, sa] = pa
        pred[b, sb] = pb

    def insere(a, b, sa, pb):
        """Insère l'arête (a, b), avec sa et pb les points tels qu'à la fin de
        l'opération, succ[a, b] = sa et pred[b, a] = pb.
        Conserve les invariants des structures de données pred et succ."""
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

    def common_tangent(x0, y0):
        """Si nous avons calculé la tangente commune à deux ensembles de points
        X et Y séparés par une droite, et que x0 et y0 sont les points de X et Y
        les plus proches de la droite, calcule la tangente commune des
        enveloppes convexes de X et Y passant par les points x de X et y de Y
        et telle que tous les points de X et de Y se trouvent à gauche du
        segment orienté (x, y)."""
        ##x = max(left)  #le point le plus à droite de l'ensemble 'left'.
        ##y = min(right) #le point le plus à gauche de l'ensemble 'right'.
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
        first[x] = y
        
        while True: 
            if orientation(x, y, pred[y, x]) == DIRECT:
                y1 = pred[y, x]
                y2 = pred[y, y1]
                while position_cercle_circonscrit(x, y, y1, y2) == DEDANS:
                    delete(y, y1)
                    y1 = y2
                    y2 = pred[y, y1]
            else:
                y1 = None


            if orientation(x, y, succ[x, y]) == DIRECT:
                x1 = succ[x, y]
                x2 = succ[x, x1]
                while position_cercle_circonscrit(x, y, x1, x2) == DEDANS:
                    delete(x, x1)
                    x1 = x2
                    x2 = succ[x, x1]
            else:
                x1 = None

            if x1 is None and y1 is None:
                break
            elif x1 is None:
                #on trace l'arête (x, y1)
                insere(y1, x, y, y)
                y = y1
            elif y1 is None:
                #on trace l'arête (y, x1)
                insere(y, x1, x, x)
                x = x1
            elif position_cercle_circonscrit(x, y, y1, x1) == DEDANS:
                insere(y, x1, x, x)
                x = x1
            else:
                insere(y1, x, y, y)
                y = y1

        #Finalement, on met à jour les informations sur l'enveloppe convexe
        #commune.
        first[y] = x
        assert (x, y) in succ
        assert (y, x) in succ

    def compute(points):
        """Calcule la triangulation de Delaunay et l'enveloppe convexe d'une
        liste de points en mettant à jour les structures de données
        succ, pred et first.

        Préconditions :
        -la liste points contient au moins deux éléments
        -les points sont tous distincts"""
        n = len(points)
        
        if n == 2:
            [a, b] = points
            #S'il n'y a que deux points, on trace le segment [a, b].
            succ[a, b] = pred[a, b] = b
            succ[b, a] = pred[b, a] = a
            first[a] = b
            first[b] = a
            
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
            elif orientation(a, b, c) == INDIRECT:
                succ[a, b] = succ[b, a] = pred[a, b] = pred[b, a] = c
                succ[a, c] = succ[c, a] = pred[a, c] = pred[c, a] = b
                succ[b, c] = succ[c, b] = pred[b, c] = pred[c, b] = a
                first[a] = c
                first[b] = a
                first[c] = b
            else:
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
                
        else: #S'il y a au moins 4 points :
            var_x, var_y = variance_xy(points)
            if var_y < var_x:
                #On sépare les points selon une droite verticale :
                med = pseudo_mediane(points)
                left = [p for p in points if p < med]
                right = [p for p in points if p >= med]
                compute(left)
                compute(right)
                x, y = common_tangent(max(left), min(right))
                merge(x, y)
            else:
                #On sépare les points selon une droite horizontale :
                med = pseudo_mediane(points, key=inv)
                down = [p for p in points if inv(p) < inv(med)]
                up = [p for p in points if inv(p) >= inv(med)]
                compute(down)
                compute(up)
                x, y = common_tangent(max(down, key=inv), min(up, key=inv))
                merge(x, y)
            
    compute(points)
    return succ
