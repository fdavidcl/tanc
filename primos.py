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
#     2. Primalidad y factorización            #primos
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
# 2. Primalidad y factorización                    #primos
##########################################################

from collections import Counter
import random
import sympy
import math

# Decorador para repetición de tests de pseudoprimalidad
def repeat_test(test):
    def repeated(self, n, k = 1):
        bases = []

        for _ in xrange(k):
            b, val = test(n)
            bases.append(b)

            if not val:
                return (bases, False)

        print("Tras {} intentos, es posible que {} sea primo".format(k, n))
        return (bases, True)
    return repeated

class PseudoprimalityTest:
    @repeat_test
    def test(self, n):
        # Descarto el caso 2
        if n == 2:
            return (2, True)
        
        # Elección de la base
        b = random.randrange(2, n)
        # Máximo común divisor
        d = sy.gcd(b, n)
        
        if d > 1:
            print("{} es divisor de {}".format(d, n))
            return (b, False)
        else:
            # Devuelve el valor de verdad de la condición
            return (b, pow(b, n-1, n) != 1)

class EulerPTest:
    def jacobisym(self, m, n):
        m = m % n
    
        if n % 2 == 0:
            raise ValueError, "n no puede ser múltiplo de 2"
    
        if gcd(m, n) != 1:
            return 0
        elif m == 1:
            return 1
        elif m % 2 == 0:
            return int((-1) ** ((n ** 2 - 1)/8) * jacobisym(m / 2, n))
        else:
            return int((-1) ** (((m - 1)*(n - 1))/4) * jacobisym(n, m))

    @repeat_test
    def epsp(self, n):
        # Se descarta el caso par
        if n % 2 == 0:
            print("{} es par".format(n))
            return (2, n == 2)
    
        # Selección de la base
        b = random.randrange(2, n)
        d = sy.gcd(b, n)
        
        if d == 1:
            # Devuelve el valor de verdad de la condición
            return (b, pow(b, int((n - 1)/2), n) == self.jacobisym(b, n) % n)
        else:
            print("{} es divisor de {}".format(d, n))
            return (b, False)

class StrongPTest:
    def mpot(self, p, n):
        return 0 if n % p != 0 else 1 + self.mpot(p, n/p)
    
    def mpot_iter(self, p, n):
        i = 0
        
        while n % p == 0:
            i += 1
            n = n / p
            
        return i

    @repeat_test
    def fpsp(self, n):
        # Descarta el caso 1
        if n == 1:
            return (1, False)
        
        # Descarta el caso par
        if n % 2 == 0:
            print("{} es par".format(n))
            return (2, n == 2)
        
        # Calcula s y t
        s = self.mpot(2, n - 1)
        t = int((n - 1) / pow(2, s))
        b = random.randrange(2, n)
        
        gcd = sy.gcd(b, n)
        if gcd != 1:
            print("{} es divisor de {}".format(gcd, n))
            return (b, False)
        
        b_t_modn = pow(b, t, n)
        if b_t_modn == 1 or b_t_modn == n - 1:
            return (b, True)
        
        for i in xrange(1, s):
            b_t_2i_modn = pow(b_t_modn, pow(2, i), n)
            if b_t_2i_modn == n - 1:
                print("{}^({} * 2^{}) = -1 mod {}".format(b, t, i, n))
                return (b, True)
            elif b_t_2i_modn == 1:
                return (b, False)
    
        return (b, False)


def compare_tests(n, k = 1):
    ps = PseudoprimalityTest().test(n, k)
    euler = EulerPTest().test(n, k)
    fuerte = StrongPTest().test(n, k)
    sym = sy.isprime(n)
    
    print("Sympy: {}\nPseudoprimalidad: {} (bases {})\nEuler-pseudoprimalidad: {} (bases {})\nPseudoprimalidad fuerte: {} (bases {})".format(
        sym,
        ps[1], ps[0],
        euler[1], euler[0],
        fuerte[1], fuerte[0]
    ))

class BaseFactor:
    def __init__(self, n):
        self.n = n
    
    def abmod(self, x):
        return x % self.n - (0 if x % self.n <= n/2 else n)

    def mpot(self, p, x):
        return 0 if x % p != 0 else 1 + self.mpot(p, x/p)
    
    def mayorpot(self, p, x):
        if p == -1:
            return int(x < 0)
        return self.mpot(p, x)

    def bnumer(self, b, base, alfav = []):
        factorizable = self.abmod(b * b)
        
        # Descarta el caso 0
        if factorizable == 0:
            return False
        
        # Divide por la mayor potencia de cada primo de la base
        for f in base:
            p = self.mayorpot(f, factorizable)
            # Acumula los exponentes en el alfa-vector
            alfav.append(p)
            factorizable = factorizable / pow(f, p)
            
        # Comprueba si hemos terminado de factorizar
        return factorizable == 1

    def vec_alfa(self, b, base):
        exp = []
        isb = self.bnumer(b, base, exp)
        return exp if isb else []

    def parp(self, l):
        return all(x % 2 == 0 for x in l)

    # ssuma_dos suma los valores de dos listas
    def ssuma_dos(self, lista1, lista2):
        if len(lista1) != len(lista2):
            raise ValueError, "Listas de distinta longitud"
        
        return map(lambda (a, b): a + b, zip(lista1, lista2))

    # ssuma suma una cantidad arbitraria de listas
    def ssuma(self, *listas):
        return reduce(self.ssuma_dos, listas)

    def selections(self, amount):
        selected = [False for _ in xrange(amount)]
        
        while any(not x for x in selected):
            i = len(selected) - 1
            
            while selected[i]:
                selected[i] = False
                i -= 1
                
            selected[i] = True
        
            # Proporciona una nueva selección. La siguiente selección no se 
            # calcula hasta que se "pida" al objeto generador
            yield selected

    def bi(self, k_max, i_max, base):
        # Lista con las raíces
        raices = [int(math.floor(math.sqrt(k * self.n))) for k in xrange(1, k_max + 1)]
        # Lista con raíces y sumas
        todos = [num + i for num in raices for i in xrange(0, i_max + 1)]
        
        # Proceso de filtrado: creo un diccionario vacío y añado los B-números
        # y sus alfavectores
        bn = {}
        
        for b in todos:
            alfa = self.vec_alfa(b, base)
            
            if alfa: # aprovecho que [] es un valor de falsedad
                bn[b] = alfa
                
        return bn

    def sol_base(self, base, bnumsdict):
        bnums = bnumsdict.keys()
        alfas = bnumsdict.values()
        
        num_select = 1
        select = self.selections(len(bnums))
        remaining = True
        
        for selected in iter(select):
            alfasum = self.ssuma(*[alfas[i] for i in xrange(len(bnums)) if selected[i]])
            
            if self.parp(alfasum):
                t = reduce(
                    lambda prod, i: (prod * (bnums[i] if selected[i] else 1)) % self.n,
                    xrange(len(bnums)),
                    1
                )
                s = reduce(
                    lambda prod, v: (prod * pow(base[v], int(alfasum[v]/2), n)) % self.n,
                    xrange(len(alfasum)),
                    1
                )
                if t != s and (-t % self.n) != s:
                    return (t, s)

    def soleq(self, primos, k_max, i_max):
        base = [-1] + [sy.generate.prime(i + 1) for i in xrange(primos)]
        bnumsdict = self.bi(k_max, i_max, base)
        return self.sol_base(base, bnumsdict)

    def soleqfc(self, s):
        frac = sy.continued_fraction_periodic(0, 1, self.n)
        l1 = [frac[0]] + frac[1] if len(frac) > 1 else frac
        s = 1 + s % len(l1)
        l2 = sy.continued_fraction_convergents(l1[:s])
    
        pbs = map(lambda x: x.p, l2)
        ab = map(lambda b: abmod(b*b, n), pbs)
        fac_ab = map(factorint, ab)
        
        factors = [f for d in fac_ab for f in d.keys()]
        exps = [e for d in fac_ab for e in d.values()]
        
        which_twice = list(set(x for x in factors if factors.count(x) > 1 and x > 0))
        which_evenexp = list(set(
            map(lambda i: factors[i], filter(lambda i: exps[i] % 2 == 0, range(len(exps))))
        ))
        
        base = [-1] + list(set(which_twice + which_evenexp))
        bnumsdict = {}
        
        for b in pbs:
            alfa = self.vec_alfa(b, base)
            if alfa:
                bnumsdict[b] = alfa
    
        return self.sol_base(base, bnumsdict)

    def fac_with(self, fn, *sol_args):
        factorizando = [self.n]
        factorizacion = Counter()
        factorizacion[self.n] += 1
        
        # Comprobaciones previas: primalidad
        if self.fpsp(self.n, 10)[1]:
            print("{} parece primo. Saliendo...".format(self.n))
            return factorizacion
        
        # Eliminamos la mayor potencia de 2
        mp2 = self.mayorpot(self.n, 2)
        quot = int(self.n/pow(2, mp2))
        if mp2 > 0:
            factorizacion[2] = mp2
            factorizacion[quot] = factorizacion[n]
            del factorizacion[n]
            
            factorizando.pop()
            factorizando.append(quot)
            
        # Proceso de factorización
        while factorizando:
            m = factorizando.pop()
            
            # Trato aparte el caso de los cuadrados
            sqrtm = int(sqrt(m))
            
            if sqrtm != 1 and m % sqrtm == 0:
                # Añade el exponente (doble del exponente de m) a la factorización de n
                factorizacion[sqrtm] = 2 * factorizacion[m]
                del factorizacion[m]
                
                # Deja sqrtm pendiente de factorizar
                factorizando.append(sqrtm)
                
                # Caso en que m no es un cuadrado
            else:
                # Cálculo de la solución (t, s)
                sol = fn(m, *sol_args)
                
                if sol:
                    t, s = sol
                    factor = sy.gcd(t - s, m)
                    mp = mayorpot(factor, m)
                    remaining = int(m/pow(factor, mp))
                    
                    # Añade el exponente del factor y del cociente de 'm' por 'factor'
                    factorizacion[factor] += factorizacion[m] * mp
                    factorizacion[remaining] += factorizacion[m]
                    del factorizacion[m]
                    
                    # Deja pendientes el factor y el cociente
                    factorizando.append(int(factor))
                    factorizando.append(int(remaining))
    
        # Eliminamos posible aparición del 1 en la factorización
        del factorizacion[1]
        return factorizacion

    def fac(self, primos, k_max = 5, i_max = 5):
        return fac_with(self.soleq, primos, k_max, i_max)
            
    def facfc(self, s = 20):
        return fac_with(self.soleqfc, s)

