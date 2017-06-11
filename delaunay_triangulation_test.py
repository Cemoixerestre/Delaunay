from delaunay_triangulation import *

def test_triangulation(points, triangulation):
    """"Teste si la triangulation de Delaunay d'un ensemble de points est
    correcte.

    Plus précisément, il s'agit ici d'une liste de points sous la forme (x, y)
    avec x et y des entiers. Cette liste contient au moins deux points, et tous
    les points de la liste sont distincts.

    On vérifie que pour tous les triangles de la triangulation, aucun point ne
    se trouve stricement à l'intérieur du cercle circonscrit au triangle."""
    for (a, b), c in triangulation.items():
        if (b, c) not in triangulation:
            continue
        #À cet endroit là, on sait que les points a, b, c sont deux à deux
        #reliés, ils forment un triangle.

        #Vérifie que le triangle a, b, c n'est pas plat :
        assert orientation(a, b, c) != ALIGNES

        d = triangulation[c, b]
        if orientation(a, b, c) == DIRECT:
            assert position_cercle_circonscrit(a, b, c, d) != DEDANS
        if orientation(a, b, c) == INDIRECT:
            assert position_cercle_circonscrit(a, c, b, d) != DEDANS

from random import randrange, sample

#« Test du test » : on vérifie que le test passe pour une triangulation
#correcte et échoue pour une triangulation incorrecte.
points = [(0, 0), (0, 1), (1, 0), (2, 1)]
[a, b, c, d] = points
#triangulation correcte :
t1 = {(a, b): c, (a, c): b, (b, a): c, (b, c): d, (b, d): a, (c, a): d,
      (c, d): b, (c, b): a, (d, b): c, (d, c): b}
#triangulation incorrecte :
t2 = {(a, b): c, (a, c): d, (a, d): b, (b, a): d, (b, d): a, (c, a): d,
      (c, d): a, (d, a): c, (d, c): b, (d, b): a}

test_triangulation(points, t1)
try:
    test_triangulation(points, t2)
    #on renvoie une erreur si aucune erreur n'a été renvoyée :
    raise Exception("Le test aurait du échouer.")
except AssertionError:
    pass

#Tests avec des points aléatoires :
def genere(n, borne, points_depart=set()):
    """Génère n points aléatoires dans le plan, chacun avec des cordonnées
    entières comprises entre 0 en borne-1"""
    points = points_depart.copy()
    x, y = randrange(borne), randrange(borne)
    for i in range(n):
        while (x, y) in points:
            x, y = randrange(borne), randrange(borne)
        points.add((x, y))
    return list(points)

for nb_points in (3, 4, 5, 10, 100, 1000, 10000, 100000):
    print("Test avec {} points aléatoires".format(nb_points))
    points = genere(nb_points, 10000)
    print("Points générés.")
    triangulation = delaunay_triangulation(points)
    print("Triangulation effectuée.")
    test_triangulation(points, triangulation)

#Tests avec des points alignés :
lx = sample(list(range(10000)), 1000)
points = [(x, 0) for x in lx]
test_triangulation(points, delaunay_triangulation(points))

lx = sample(list(range(10000)), 1000)
points = [(x, x) for x in lx]
test_triangulation(points, delaunay_triangulation(points))

#Tests avec des points alignés et des points non alignés :
points0 = set(points)
for n in (100, 1000, 10000):
    points = genere(n, 10000, points0)
    test_triangulation(points, delaunay_triangulation(points))
    
print("Tous les tests ont été passés avec succès.")
