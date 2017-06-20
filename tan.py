#!/usr/bin/env python2
# -*- coding: utf-8 -*-

from __future__ import print_function, division, generators, absolute_import
import sympy as sy
from collections import Counter
from PIL import Image

##########################################################
# Teoría de Numeros y Criptografía
# Prácticas
#
# Francisco David Charte Luque
#---------------------------------------------------------
#
# Índice
#   TEORÍA DE NÚMEROS                          #tanc
#     3. Factorización en dominios euclídeos   #facde
#     4. Factorizaión de ideales               #facideal
#     5. Cálculo del grupo de clase            #classgr
#
# También se puede buscar la implementación de cualquier
# método de una clase como Clase#metodo, por ejemplo
# Ideal#norma.
# 
##########################################################

##########################################################
# TEORÍA DE NÚMEROS                                   #tan
##########################################################    

##########################################################
# 3. Factorización en dominios euclídeos            #facde
##########################################################
    
def libre_de_cuadrados(n):
    return all(e == 1 for e in sy.ntheory.factor_.factorint(n).values())

class Numero:
    """
    Clase que envuelve un número de un cuerpo de números Q(√d). Es conveniente utilizar el 
    método `get` de la clase AnilloEnteros para obtener instancias de esta clase, en lugar
    de crearlas directamente
    """
    def __init__(self, num, anillo):
        """
        Numero#Numero(num : sy.Number, anillo : AnilloEnteros) : Numero
        Constructor
        """
        self.num = sy.sympify(num).expand()
        self.anillo = anillo

        # Precálculo de algunas constantes
        # self._x = self.num.coeff(self.anillo.sqrtd, 0)
        self._y = self.num.coeff(self.anillo.sqrtd)
        self._x = (self.num - self._y * self.anillo.sqrtd).expand()
        self._a = (self._x - self._y) if self.anillo.es1mod4 else self._x
        self._b = (2 * self._y) if self.anillo.es1mod4 else self._y
        self._conjugado = self._x - self._y * self.anillo.sqrtd
        self._norma = (self._x * self._x - self.anillo.d * self._y * self._y)
        self._traza = 2 * self._x
        self._es_entero = self._norma.is_Integer and self._traza.is_Integer
        
    def __str__(self):
        """
        Numero#str() : string
        Traduce el número a string
        """
        return "{}{}{}{}".format(
            self.x() if self.x() != 0 else "",
            "+" if self.y() > 0 and self.x() != 0 else "",
            self.y() if self.y() != 0 else "",
            self.anillo.str_d if self.y() != 0 else ""
        )

    def describe(self):
        """
        Numero#describe() : void
        Imprime por pantalla una descripción del número
        """
        print("{} en Q({}), {}es entero de {}".format(
            str(self.num), self.anillo.str_d, "" if self.es_entero() else "no ",
            str(self.anillo)))
    
    def norma(self):
        """
        Numero#norma() : sy.Number
        Calcula la norma
        """
        return self._norma

    def traza(self):
        """
        Numero#traza() : sy.Number
        Calcula la traza
        """
        return self._traza

    def es_entero(self):
        """
        Numero#es_entero() : bool
        Indica si el número es un entero del anillo
        """
        return self._es_entero

    def x(self):
        """
        Numero#x() : sy.Number
        Calcula la componente x en x + y sqrt(d)
        """
        return self._x

    def y(self):
        """
        Numero#y() : sy.Number
        Calcula la componente y en x + y sqrt(d)
        """
        return self._y
    
    def xy(self):
        """
        Numero#xy() : (sy.Number, sy.Number)
        Calcula x e y en x + y sqrt(d)
        """
        return self.x(), self.y()

    def a(self):
        """
        Numero#a() : sy.Number
        Calcula la componente a en a + b e
        """
        return self._a

    def b(self):
        """
        Numero#b() : sy.Number
        Calcula la componente b en a + b e
        """
        return self._b
    
    def ab(self):
        """
        Numero#ab() : (sy.Number, sy.Number)
        Calcula a y b en a + b e
        """
        return self.a(), self.b()

    def conjugado(self):
        """
        Numero#conjugado() : Numero
        Calcula el conjugado y lo devuelve
        """
        return self.anillo.get(self._conjugado)

    def divide(self, b):
        """
        Numero#divide(b : Numero) : bool
        Indica si el número divide a b
        """
        if b.anillo.d != self.anillo.d:
            raise ValueError
        return b.cociente(self).es_entero()

    def cociente(self, b):
        """
        Numero#cociente(b : Numero) : Numero
        Calcula el cociente self / b
        """
        if b.anillo.d != self.anillo.d:
            raise ValueError
        return Numero(self.num * b.conjugado().num / b.norma(), self.anillo)

    def es_unidad(self):
        """
        Numero#es_unidad() : bool
        Indica si el número es una unidad
        """
        return self.es_entero() and self.norma() == 1

    def es_irreducible(self):
        """
        Numero#es_irreducible() : bool
        Indica si el número es un irreducible
        """
        nor = self.norma()
        nsq = sy.sqrt(abs(nor))
        return nor.is_prime or (nsq.is_prime and len(list(self.anillo.connorma(nsq))) == 0)

    def es_asociado(self, b):
        """
        Numero#es_asociado(b : Numero) : bool
        Indica si el número es un asociado de b
        """
        return self.divide(b) and b.divide(self)

    def factoriza(self):
        """
        Numero#factoriza() : Dict[Numero, int]
        Factoriza el número y devuelve un diccionario donde las claves son 
        los factores y los valores son los exponentes
        """
        result = Counter()
        alpha = self

        # memoizamos la búsqueda de elementos por norma, que puede
        # ser lenta en el caso d > 0
        cn = dict()

        # En la práctica 4 faltaba la comprobación de que alpha
        # no fuese unidad, lo que daba lugar a problemas (e.g.
        # factorizar 3 en O(Q(sqrt(3))) )
        while not alpha.es_irreducible() and not alpha.es_unidad():
            fact = sy.ntheory.factor_.factorint(alpha.norma())
            
            p = fact.keys()[0]

            if not (p in cn.keys()):
                cn[p] = list(self.anillo.connorma(p))

                if self.anillo.d > 0:
                    cn[p] += list(self.anillo.connorma(-p))
            
                # Si no hay elementos de norma p usaremos el propio p
                cn[p].append(self.anillo.get(p))

            for a in cn[p]:
                c = alpha.cociente(a)
                if c.es_entero():
                    result[a.num] += 1
                    alpha = c
                    break
        
        if not alpha.es_unidad():
            result[alpha.num] += 1
    
        return result

class AnilloEnteros:
    def __init__(self, d):
        """
        AnilloEnteros#AnilloEnteros(d : int) : AnilloEnteros
        Constructor del anillo de enteros. Decide si d = 1 mod 4 y demás
        """
        self.d = d
        self.sqrtd = sy.sqrt(self.d)
        self.str_d = "√{}".format(self.d)

        self.es_dominio_euclideo = d in [-1, -2, -3, -7, -11,
                                         2, 3, 5, 6, 7, 11, 13,
                                         17, 19, 21, 29, 33, 37,
                                         41, 57, 73]
        
        self.es1mod4 = (d % 4) == 1

        if self.es1mod4:
            self.e = sy.Rational(1, 2) + sy.Rational(1, 2) * sy.sqrt(d)
            self.str_e = "(1 + {}) / 2".format(self.str_d)
        else:
            self.e = sy.sqrt(d)
            self.str_e = self.str_d

        # Memoizamos los números: intentamos mantener una sola copia
        # de cada número, ahorrando así muchos cálculos innecesarios
        self.numeros = dict()
        # Precalculamos las soluciones de Pell
        self._pell1_solutions = self.pell() if d > 0 else []

    def __str__(self):
        """
        AnilloEnteros#str() : string
        Da la string que identifica al anillo
        """
        return "O(Q({})) = Z({})".format(self.str_d, self.str_e)
    
    def describe(self):
        """
        AnilloEnteros#describe()
        Imprime por pantalla una descripción del anillo
        """
        print("Anillo de enteros de Q({}), Z({})".format(self.str_d, self.str_e))

    def __eq__(self, otro):
        """
        AnilloEnteros#==(otro : AnilloEnteros) : bool
        Indica si dos anillos de enteros son iguales
        """
        return self.d == otro.d
        
    def get(self, num):
        """
        AnilloEnteros#get(num : sy.Number) : Numero
        Obtiene un número del cuerpo de números (posiblemente del anillo)
        """
        if not num in self.numeros:
            self.numeros[num] = Numero(num, self)
        
        return self.numeros[num]

    def fast_is_square(self, positive):
        """
        AnilloEnteros#fast_is_square(positive : int) : bool

        Indica si un entero es un cuadrado perfecto. Mucho más rápido que calcular la raíz cuadrada 
        simbólica en sympy y comprobar si esa es entera.
        Fuente: https://stackoverflow.com/a/2489519
        """
        if positive <= 3:
            return positive == 0 or positive == 1
        x = positive // 2
        seen = set([x])
        while x * x != positive:
            x = (x + (positive // x)) // 2
            if x in seen:
                del(seen)
                return False
            seen.add(x)
        del(seen)
        return True
        
    
    def pell_solutions_for_limits(self, n, inf, sup):
        """
        AnilloEnteros#pell_solutions_for_limits(n : int, inf : float, sup : float) : Generator[(int, int)]
        Lista las soluciones de la ecuación de pell entre los límites indicados
        """
        sols = set()

        if sup - inf > 1e4: # avisamos cuando vayamos a hacer más de 10000 iteraciones
            print("[Aviso Pell: {} iteraciones]".format(sup - inf))
    
        for y in xrange(inf, sup + 1):
            x = int(self.d * y * y + n)
            
            if self.fast_is_square(x):
                x = sy.sqrt(x)
                for s in [(x, y), (-x, y), (x, -y), (-x, -y)]:
                    if not s in sols:
                        sols.add(s)
                        yield s
        del(sols)
    
    def eqpell_dneg(self, n):
        """
        AnilloEnteros#eqpell_dneg(n : int) : Generator[(int, int)]
        Encuentra las soluciones de la ecuación de Pell para d < 0
        """
        if n < 0:
            return set()
        limit = sy.sqrt(sy.Rational(n, -self.d))
        return self.pell_solutions_for_limits(n, 0, limit)

    def pell(self):
        """
        AnilloEnteros#pell_dneg() : List[(int, int)]
        Encuentra las soluciones de la ecuación de Pell para d > 0 y n = 1
        """
        contf = sy.continued_fraction_periodic(0, 1, self.d)

        l = []
        if len(contf) > 0:
            l.append(contf[0])
            if len(contf) > 1:
                l += contf[1][:-1]
                
        convg = list(sy.continued_fraction_convergents(l))[-1]
        x0, y0 = convg.p, convg.q
        
        if x0 * x0 - self.d * y0 * y0 == 1:
            return (x0, y0)
        else:
            return (x0 * x0 + y0 * y0 * self.d, 2 * x0 * y0)

    def eqpell_dpos(self, n):
        """
        AnilloEnteros#eqpell_dpos(n : int) : Generator[(int, int)]
        Encuentra las soluciones de la ecuación de Pell para d > 0
        """
        r, s = self._pell1_solutions
            
        inf, sup = None, None

        if n == 1:
            return set([(r, s)])
        elif n > 0:
            inf = 0
            sup = int(sy.sqrt(sy.Rational(n * (r - 1), (2 * self.d))))
        else:
            inf = int(sy.sqrt(sy.Rational(-n, self.d)))
            sup = int(sy.sqrt(sy.Rational(-n * (r + 1), (2 * self.d))))

        return self.pell_solutions_for_limits(n, inf, sup)

    def eqpell(self, n):
        """
        AnilloEnteros#eqpell(n : int) : Generator[(int, int)]
        Encuentra las soluciones de la ecuación de Pell para cualquier d
        """
        return self.eqpell_dpos(n) if self.d > 0 else self.eqpell_dneg(n)

    def connorma_xmod4(self, n):
        """
        AnilloEnteros#connorma_xmod4(n : int) : Generator[Numero]
        Encuentra elementos con norma n para d != 1 mod 4
        """
        for (x, y) in iter(self.eqpell(n)):
            yield self.get(x + y * sy.sqrt(self.d))

    def connorma_1mod4(self, n):
        """
        AnilloEnteros#connorma_1mod4(n : int) : Generator[Numero]
        Encuentra elementos con norma n para d == 1 mod 4
        """
        for (x, y) in iter(self.eqpell(4 * n)):
            if (x - y) % 2 == 0:
                yield self.get((x + (y * sy.sqrt(self.d))) / 2)
    
    def connorma(self, n):
        """
        AnilloEnteros#connorma(n : int) : List[Numero]
        Encuentra elementos con norma n
        """
        if self.es1mod4:
            return self.connorma_1mod4(n)
        else:
            return self.connorma_xmod4(n)

    def irreducible_e(self):
        """
        AnilloEnteros#irreducible_e() : 3-tuple[int]

        Devuelve los coeficientes del polinomio irreducible asociado a e
        """
        # nor(e) - tr(e) X + X²
        return (self.get(self.e).norma(), -self.get(self.e).traza(), 1)
    
    def ramifica(self, p):
        """
        AnilloEnteros#ramifica(p : int) : bool

        Indica si el primo ramifica en el anillo de enteros
        """
        return len(self.raices_e(p)) > 1

    def raices_e(self, p):
        """
        AnilloEnteros#raices_e(p : int) : List[int]

        Encuentra las raíces del irreducible asociado a e en módulo p
        """
        if not sy.sympify(p).is_prime:
            raise ValueError, "Necesito un primo!"

        # p ramifica <=> Irr(e)_p tiene raíces en Z_p
        irr = self.irreducible_e()

        raices = []
        # recorro Z_p buscando una raíz
        for i in xrange(p):
            if ((irr[0] + irr[1] * i + irr[2] * i * i) % p) == 0:
                raices.append(i)
        
        if len(raices) == 1:
            raices.append(raices[0])
            
        if not raices:
            raices = [p]

        return raices

    def ideal_cero(self):
        """
        AnilloEnteros#ideal_cero() : Ideal

        Construye el ideal cero del anillo
        """
        return Ideal(self.get(0))

    def ideal_total(self):
        """
        AnilloEnteros#ideal_total() : Ideal

        Ve el anillo como un ideal de sí mismo
        """
        return Ideal(self.get(1))

    def ideales_divisores(self, p):
        """
        AnilloEnteros#ideales_divisores(p : int) : List[Ideal]
        Calcula los ideales divisores de un primo p
        """
        divisores = []
        
        if self.ramifica(p):
            x1, x2 = self.raices_e(p)
            divisores.append(
                Ideal([self.get(p), self.get(self.e - x1)])
            )
            divisores.append(
                Ideal([self.get(p), self.get(self.e - x2)])
            )
        else:
            divisores.append(Ideal(self.get(p)))

        return divisores

    def cota_minkowski(self):
        """
        AnilloEnteros#cota_minkowski() : sympy.Number
        Proporciona la cota de Minkowski para este anillo
        """

        s = 2 if self.d > 0 else 0
        t = 1 if self.d < 0 else 0
        n = s + 2 * t
        disc = self.d * (1 if self.es1mod4 else 4)
        mst = (4 / sy.pi) ** t * (sy.factorial(n) / (n ** n))
        return mst * sy.sqrt(abs(disc))


##########################################################
# 4. Factorización de ideales                    #facideal
##########################################################

"""
Claves para usar en métodos de ordenación
"""
def _cmp_key0(par):
    return float("inf") if par[0] == 0 else abs(par[0])

def _cmp_key1(par):
    return float("inf") if par[1] == 0 else abs(par[1])

class Ideal:
    """
    Clase que representa un ideal de un anillo de enteros. Permite 
    factorizar ideales de anillos de enteros que no sean dominios
    euclídeos, entre otras funcionalidades
    """
    
    def __init__(self, generadores):
        """
        Ideal#Ideal(generadores : Numero | List[Numero]) : Ideal

        Constructor de Ideales
        """
        if not generadores:
            raise ValueError, "Necesito al menos un generador (el Numero 0 para el ideal cero)"

        if isinstance(generadores, Numero):
            generadores = [generadores]
        
        self.anillo = generadores[0].anillo
        self.generadores = self.generadores_iniciales = generadores

        self.mat = self.matriz_relatores()
        self.relatores = self.reduce_relatores()
        
        a, b = self.relatores[0]
        _, c = self.relatores[1]
            
        # Guardamos ambos generadores, sabemos que son del tipo g1 = a + be y g2 = ce
        g1 = self.anillo.get(a + self.anillo.e * b)
        g2 = self.anillo.get(self.anillo.e * c)
        self.generadores = [g1, g2]

    def __str__(self):
        """
        Ideal#str() : string
        """
        return "<{}, {}>".format(str(self.generadores[0]), str(self.generadores[1]))

    def describe(self):
        """
        Ideal#describe() : void

        Imprime por pantalla una descripción del ideal
        """
        print("Ideal generado por: {}, {}. Norma: {}{}\nMatriz de relatores reducida:\n/{} {}\\\n\\{} {}/\n".format(
            str(self.generadores[0]), str(self.generadores[1]), self.norma(),
            ", ideal primo" if self.es_primo() else "",
            self.relatores[0][0], self.relatores[1][0],
            self.relatores[0][1], self.relatores[1][1]
        ))
        
    def matriz_relatores(self):
        """
        Ideal#matriz_relatores() : List[List[int]]

        Calcula la matriz de relatores del grupo (I, +) a partir de los generadores del ideal I
        """
        generadores_grupo = self.generadores_iniciales + map(
            lambda x: x.anillo.get(sy.expand(x.anillo.e * x.num)),
            self.generadores_iniciales)
        return map(lambda x: list(x.ab()), generadores_grupo)

    def lattice_reduce(self):
        """
        Ideal#lattice_reduce() : void

        Reduce la matriz de relatores mediante operaciones elementales
        de columnas
        """

        while True:
            if any(col[0] != 0 for col in self.mat[1:]):
                #self.mat = sorted(self.mat, key = self._cmp_key(0))
                self.mat.sort(key = _cmp_key0)
            
                if self.mat[0][0] == 0:
                    raise ValueError, "No puedo reducir esta matriz"

                for col in xrange(1, len(self.mat)):
                    q = self.mat[col][0] // self.mat[0][0]
                    self.mat[col][0] -= q * self.mat[0][0]
                    self.mat[col][1] -= q * self.mat[0][1]
                
            elif any(col[1] != 0 for col in self.mat[2:]):
                self.mat[1:] = sorted(self.mat[1:], key = _cmp_key1)
            
                if self.mat[1][1] == 0:
                    raise ValueError, "No puedo reducir esta matriz"

                # En este momento todos los valores de la primera fila son 0 (salvo mat[0][0])
                for col in xrange(2, len(self.mat)):
                    self.mat[col][1] = self.mat[col][1] % self.mat[1][1]
                
            elif self.mat[1][1] != 0:
                # reduzco b dejando el resto por c
                self.mat[0][1] = self.mat[0][1] % self.mat[1][1]

                # Termina la ejecución
                return
            else:
                return
                

    def reduce_relatores(self):
        """
        Ideal#generadores_minimos() : List[List[int]]

        Devuelve la matriz 2x2 de generadores mínimos a partir de la 
        reducción de la matriz de relatores
        """
        self.lattice_reduce()

        if any((e[0] != 0 or e[1] != 0) for e in self.mat[2:]):
            raise AssertionError, "Hay ceros donde no deberia"
        
        return self.mat[:2]

    def norma(self):
        """
        Ideal#norma() : sy.Number
        
        Calcula la norma del ideal
        """
        return abs(sy.sympify(self.relatores[0][0] * self.relatores[1][1]))

    def es_elemento(self, numero):
        """
        Ideal#es_elemento(numero : Numero) : bool

        Decide si un número del anillo pertenece al ideal
        """
        assert(numero.anillo == self.anillo)

        a, b = self.relatores[0]
        _, c = self.relatores[1]
        n, m = numero.ab()

        # Caso en que no podemos intentar dividir por a y c
        # a == 0 <=> c == 0 <=> I = 0
        if self.norma() == 0:
            return n == 0 and m == 0

        # Comprueba que hay soluciones enteras a las ecuaciones
        return (n % a == 0) and ((m - (b * n) / a) % c == 0)

    def contiene(self, otro):
        """
        Ideal#contiene(otro : Ideal) : bool
        
        Decide si otro está contenido en el ideal
        """
        assert(otro.anillo == self.anillo)

        # Comprueba si cada generador pertenece
        return all(self.es_elemento(g) for g in otro.generadores)

    divide = contiene
    """
    Ideal#divide(otro : Ideal) : bool

    Decide si el ideal divide a otro
    """
    
    def cociente(self, otro):
        """
        Ideal#cociente(otro : Ideal) : Ideal

        Devuelve el ideal cociente
        """
        if self.anillo != otro.anillo:
            raise ValueError, "Necesito ideales del mismo anillo"
        
        if not otro.es_primo():
            raise NotImplementedError, "No puedo calcular el cociente por un ideal no primo"

        producto = self
        p = sy.sympify(otro.norma())
        
        if not p.is_prime:
            # pO es un ideal con único generador p
            p = sy.sqrt(otro.norma())
        else:
            # el divisor es un ideal primo del tipo <p, e - a>
            x1, x2 = self.anillo.raices_e(p)
            ideales = {
                x1: Ideal([self.anillo.get(p), self.anillo.get(self.anillo.e - x1)]),
                x2: Ideal([self.anillo.get(p), self.anillo.get(self.anillo.e - x2)])
            }
            
            y = x1 if otro == ideales[x2] else x2
            producto = self * ideales[y]

        return Ideal(map(lambda x: self.anillo.get(x.num / p), producto.generadores))
            
    # Operador "/" hace el cociente
    __div__ = cociente

    def producto(self, otro):
        """
        Ideal#producto(otro : Ideal) : Ideal

        Devuelve el ideal producto
        """
        if self.anillo != otro.anillo:
            raise ValueError, "Necesito ideales del mismo anillo"

        nuevos_gens = []

        for g in self.generadores:
            for h in otro.generadores:
                nuevos_gens.append(self.anillo.get(g.num * h.num))

        return Ideal(nuevos_gens)

    # Operador "*" hace el producto de ideales
    __mul__ = producto

    def producto_l(self, otros):
        """
        Ideal#producto_l(otros : List[Ideal]) : Ideal

        Versión más rápida de una multiplicación de varios ideales
        """
        nuevos_gens = self.generadores

        for i in otros:
            for g in nuevos_gens:
                more = []
                for h in i.generadores:
                    more.append(self.anillo.get(g.num * h.num))
                nuevos_gens += more

        return Ideal(nuevos_gens)

    def __pow__(self, potencia):
        """
        Ideal#**(potencia : int) : Ideal
        Calcula la potencia indicada del ideal
        """
        if potencia == 0:
            return self.anillo.ideal_total()
        return reduce(lambda ideal, _: ideal * self, xrange(1, potencia), self)

    def __eq__(self, otro):
        """
        Ideal#==(otro : Ideal) : bool

        Decide si dos ideales son iguales
        """
        if isinstance(otro, AnilloEnteros):
            # Si comparamos con el anillo, delegamos en el ideal total
            return self == otro.ideal_total()
        elif isinstance(otro, self.__class__):
            # Comprobamos igualdad por doble inclusión
            return self.contiene(otro) and otro.contiene(self)
        elif otro == 0:
            # Si comparamos con el cero lo hacemos con el ideal nulo
            return self == self.anillo.ideal_cero()
        else:
            # No son comparables
            return False

    def es_primo(self):
        """
        Ideal#es_primo() : bool
        
        Indica si el ideal es primo (equiv. maximal)
        """
        norma = self.norma()
        nsq = sy.sqrt(norma)

        # bien la norma es un primo, o bien es el cuadrado de un primo
        # que no ramifica en el anillo de enteros
        return norma.is_prime or (nsq.is_prime and not self.anillo.ramifica(nsq))
        
    def factoriza(self, debug = False):
        """
        Ideal#factoriza() : List[Ideal]

        Factoriza el ideal como producto de ideales primos
        """
        actual = self
        factorizacion = []
        while not (actual.es_primo() or actual == 0 or actual == actual.anillo):
            factnor = sy.ntheory.factor_.factorint(actual.norma())
            if debug:
                print("Norma del ideal actual {}".format(actual.norma()))
            
            p = factnor.keys()[0]
            if debug:
                print("Intentamos dividir por un ideal primo divisor de {}".format(p))
            
            for candidato in actual.anillo.ideales_divisores(p):
                if candidato.divide(actual):
                    if debug:
                        print("Dividimos {} por {}".format(actual, candidato))
                    
                    factorizacion.append(candidato)
                    actual = actual.cociente(candidato)
                    break

        if actual.es_primo():
            factorizacion.append(actual)
        return factorizacion

    def es_principal(self):
        """
        Ideal#es_principal() : bool
        Indica si el ideal está generado por un solo elemento del anillo
        """
        n = self.norma()
        l = self.anillo.connorma(n)
        for x in l:
            if self.es_elemento(x):
                return True
        l =  self.anillo.connorma(-n)
        for x in l:
            if self.es_elemento(x):
                return True
        return False
    
    def clase(self):
        """
        Ideal#clase() : Clase
        Devuelve la clase que tiene al ideal como representante
        """
        return Clase(self)

def producto_ideales(lista):
    """
    producto_ideales(lista : List[Ideal]) : Ideal
    Realiza el producto de varios ideales
    """
    return reduce(lambda i1, i2: i1 * i2, lista)

##########################################################
# 5. Cálculo del grupo de clase                   #classgr
##########################################################

class Elemento:
    """
    Representación simbólica de un elemento del grupo de clase
    """
    
    def __init__(self, bases, exponentes, ordenes):
        """
        Elemento#Elemento(bases : Ideal | List[Ideal], exponentes : int | List[int], ordenes : int | List[int]) : Elemento

        No almacena el representante de la clase, sino que contiene el producto
        de potencias de generadores que lo generan. Así, al calcular potencias y
        productos con otros elementos obtenemos siempre representantes lo más
        pequeños posibles, usando propiedades de grupos finitos
        """
        self.bases = bases if isinstance(bases, list) else [bases]
        self.exponentes = exponentes if isinstance(exponentes, list) else [exponentes]
        # Órdenes de las bases, no de las bases elevadas a los exponentes
        self.ordenes = ordenes if isinstance(ordenes, list) else [ordenes]
        
        for i in xrange(len(self.exponentes)):
            self.exponentes[i] = self.exponentes[i] % self.ordenes[i]

        self._orden = None

    def __str__(self):
        """
        Elemento#str() : string
        """
        return "".join(map(
            lambda i: str(self.bases[i]) + "^" + str(self.exponentes[i]),
            xrange(len(self.bases))
        ))

    def describe(self):
        """
        Elemento#describe() : void
        """
        print(str(self))
        
    def orden(self):
        """
        Elemento#orden() : int

        Devuelve el orden, calculándolo una vez si es necesario
        """
        if not self._orden:
            self._orden = self.calcula_orden()
        return self._orden

    def __eq__(self, otra):
        """
        Elemento#==(otra : Elemento | int) : bool
        Comprueba si dos clases son iguales. Se puede comparar una clase con el 1 para ver si es el uno del grupo
        """
        if otra == 1:
            return self.to_ideal().es_principal()
        else:
            # Buscamos el orden más bajo, para ahorrar cálculos cuando sea posible
            o1 = self.orden()
            o2 = otra.orden()
            if o1 < o2:
                return (self ** (o1 - 1)).to_ideal() * otra.to_ideal() == 1
            else:
                return (otra ** (o2 - 1)).to_ideal() * self.to_ideal() == 1

    def __pow__(self, potencia):
        """
        Elemento#**(potencia : int) : Elemento
        Calcula la potencia indicada de una clase
        """
        nu_exp = map(
            lambda e: e * potencia, self.exponentes
        )
        return Elemento(self.bases, nu_exp, self.ordenes)

    def __mul__(self, otra):
        """
        Elemento#*(otra : Elemento) : Elemento
        Multiplica dos clases simbólicamente y devuelve el producto
        """

        # no es óptimo: debería comprobar si algunas parejas de generadores son el mismo
        return Elemento(
            self.bases + otra.bases,
            self.exponentes + otra.exponentes,
            self.ordenes + otra.ordenes
        )
        
    def to_ideal(self):
        """
        Elemento#to_ideal() : Ideal
        Devuelve un ideal representante de la clase
        """
        ideal = self.bases[0] ** self.exponentes[0]
        for i in xrange(1, len(self.exponentes)):
            ideal = ideal * self.bases[i] ** self.exponentes[i]
        return ideal
    
    def calcula_orden(self):
        """
        Elemento#calcula_orden() : int
        Calcula el orden de la clase, es decir, a qué potencia se debe elevar un
        representante para que el resultante sea principal
        """
        ordenes_potencias = map(
            lambda i: self.ordenes[i] / sy.gcd(self.ordenes[i], self.exponentes[i]),
            xrange(len(self.ordenes))
        )
        if len(ordenes_potencias) == 1:
            return ordenes_potencias[0]

        limite = sy.lcm(ordenes_potencias)
        divisores = sy.divisors(limite)

        for d in divisores:
            if (self ** d).to_ideal().es_principal():
                return d

        # No deberíamos llegar aquí
        raise RuntimeError, "Nos hemos pasado de posibles órdenes"

import itertools
    
class GrupoClase:
    def __init__(self, anillo, debug = False):
        """
        GrupoClase#GrupoClase(anillo : AnilloEnteros) : GrupoClase
        """
        self.anillo = anillo
        self.debug = debug

        if not libre_de_cuadrados(self.anillo.d):
            raise ValueError, "d = {} no es libre de cuadrados!".format(self.anillo.d)
        elif debug:
            print("d = {} es libre de cuadrados".format(self.anillo.d))
        
        if debug:
            print("Calculamos las clases de ideales generadoras del grupo de clase:")
        
        self.ideales = self.calcula_generadores(debug)
                    
        self.ordenes = map(self.orden_clase, self.ideales)
        self.generadores = map(
            lambda i: Elemento(self.ideales[i], 1, self.ordenes[i]),
            xrange(len(self.ideales))
        )
        self.relaciones = set([
            tuple((j == i) * self.ordenes[i] for j in xrange(len(self.ideales)))
            for i in xrange(len(self.ideales))
        ])

    def calcula_generadores(self, debug = False):
        """
        GrupoClase#calcula_generadores() : List[Ideal]
        """
        # Calculamos la cota de Minkowski
        mink = self.anillo.cota_minkowski()
        if debug:
            print("Cota de Minkowski {} ~ {}".format(mink, float(mink)))
        # Hallamos los primos menores y nos quedamos con los que ramifican
        # Para calcular el grupo de clase basta con coger los primos por debajo de la cota de Minkowski
        # Fuente: Stein, W. (2012). Algebraic number theory, a computational approach. Harvard, Massachusetts.
        primos_menores = filter(lambda p: self.anillo.ramifica(p), sy.sieve.primerange(0, mink))
        if debug:
            print("Primos menores que la cota de Minkowski que ramifican:\n{}".format(str(primos_menores)))
        # Construimos los ideales cuyas clases serán los generadores del grupo
        generadores = map(lambda p: self.anillo.ideales_divisores(p)[0], primos_menores)
        if debug:
            if generadores:
                print("Ideales divisores de los primos:")
                for i in generadores:
                    print("\t{}".format(str(i)))
        # Nos quedamos con los generadores que no son 1
        generadores = filter(lambda g: not g.es_principal(), generadores)
        if debug:
            if generadores:
                print("Ideales divisores de los primos, no principales:")
                for i in generadores:
                    print("\t{}".format(str(i)))
            else:
                print("No se han encontrado generadores. El grupo de clase es trivial")
        return generadores

    def orden_clase(self, ideal):
        """
        GrupoClase#orden_clase(ideal : Ideal) : int
        """
        ord = 1
        pot = ideal

        while not pot.es_principal():
            ord += 1
            pot = pot * ideal

        return ord
            
    def combinations(self):
        """
        GrupoClase#combinations()
        Calcula todas las posibles relaciones entre los generadores
        """
        exponents = []
        for g in self.generadores:
            exponents.append(xrange(g.orden()))

        return itertools.product(*exponents)

    def naive_bruteforce(self, debug = True):
        """
        GrupoClase#naive_bruteforce()
        Comprueba todas las posibles relaciones entre los generadores y añade
        las que dan 0
        """
        combinations = self.combinations()
        
        for exps in iter(combinations):
            if debug:
                print("Comprobando la relación {}".format(exps))
            elemento = self.generadores[0] ** exps[0]
            
            for i in xrange(1, len(exps)):
                elemento = elemento * self.generadores[i] ** exps[i]

            if elemento == 1:
                if debug:
                    print("Es una relación del grupo!")
                self.relaciones.add(exps)

        return self.relaciones

    def pairs(self, debug = True):
        """
        GrupoClase#pairs()
        Comprueba todas las relaciones entre parejas de generadores
        """
        parejas = itertools.combinations(xrange(len(self.generadores)), 2)
        for i, j in parejas:
            g, h = self.generadores[i], self.generadores[j]
            exps = itertools.product(xrange(g.orden()), xrange(h.orden()))

            for e, f in exps:
                if g ** e * h ** f == 1:
                    if debug:
                        print("Encontrada la relación g{}^{} g{}^{} = 1".format(i, e, j, f))
                    self.relaciones.add(
                        tuple((k == i) * e + (k == j) * f for k in xrange(len(self.generadores)))
                    )

    def elimina_generador(self, i):
        """
        GrupoClase#elimina_generador(i : int) : void
        Elimina un generador y las relaciones en que aparece
        """
        del self.ideales[i]
        del self.ordenes[i]
        del self.generadores[i]
        symdiff = set()
        for r in self.relaciones:
            if r[i] != 0:
                symdiff.add(r)
            else:
                rr = list(r)
                del rr[i]
                symdiff.add(r)
                symdiff.add(tuple(rr))

        self.relaciones = self.relaciones.symmetric_difference(symdiff)
                
    def elimina_generadores_emparejados(self, debug = True):
        """
        GrupoClase#elimina_generadores_emparejados() : void
        Elimina los generadores sobrantes utilizando relaciones entre parejas 
        de generadores

        Criterios para descartar un generador:
        - si uno es potencia de otro, incluyendo los casos:
          - si son inversos
          - si son iguales
        """

        indexes = xrange(len(self.ordenes))
        a_eliminar = set()

        for r in self.relaciones:
            # Construimos la lista
            # [(g, e, o), (h, f, p)]
            # donde g, h son los índices de los generadores relacionados,
            # e, f los exponentes y o, p los órdenes.
            geo = filter(
                lambda (gen, exp, orden): exp != 0,
                zip(indexes, r, self.ordenes)
            )

            # Sólo tratamos con relaciones de pares de generadores
            if len(geo) != 2:
                continue

            g, e, o = geo[0]
            h, f, p = geo[1]

            if f == 1 and not g in a_eliminar: # h es potencia de g
                if debug:
                    print("g{} es potencia de g{}, lo eliminaremos".format(h, g))
                a_eliminar.add(h)
            elif e == 1 and not h in a_eliminar: # g es potencia de h
                if debug:
                    print("g{} es potencia de g{}, lo eliminaremos".format(g, h))
                a_eliminar.add(g)

        a_eliminar = list(a_eliminar)
                
        for i in xrange(len(a_eliminar)):
            self.elimina_generador(a_eliminar[i])

            # decrementamos los índices que apuntan a otros generadores
            for j in xrange(i, len(a_eliminar)):
                if a_eliminar[j] > a_eliminar[i]:
                    a_eliminar[j] -= 1
