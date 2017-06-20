#!/usr/bin/env python2
# -*- coding: utf-8 -*-

#  This Source Code Form is subject to the terms of the Mozilla Public
#  License, v. 2.0. If a copy of the MPL was not distributed with this
#  file, You can obtain one at http://mozilla.org/MPL/2.0/.

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
#   CRIPTOGRAFÍA                               #c
#     1. Criptología                           #cripto
#
# También se puede buscar la implementación de cualquier
# método de una clase como Clase#metodo, por ejemplo
# Ideal#norma.
# 
##########################################################

##########################################################
# CRIPTOGRAFÍA                                          #c
##########################################################

class ZnModulo:
    """
    Esta clase representa un módulo de Zn para algún valor de n.
    """
    
    def __init__(self, mod, rank):
        """
        ZnModulo#ZnModulo(mod : int, rank : int) : ZnModulo
        """
        self.mod = mod
        self.rank = rank
        
    def element(self, *scalars):
        """
        ZnModulo#element(*scalars : List[int]) : (*int)
        Devuelve un elemento del Zn-módulo a partir de los escalares de Z dados
        """
        return (scalars[i] % self.mod for i in range(self.rank))

    def sum(self, x, y):
        """
        ZnModulo#sum(x : (*int), y : (*int)) : (*int)
        Suma dos tuplas del Zn-módulo
        """
        return ((a + b) % self.mod for a, b in zip(x, y))

    def mult(self, a, x):
        """
        ZnModulo#mult(a : int, x : (*int)) : (*int)
        Multiplica un escalar por un elemento del Zn-módulo
        """
        return ((a * b) % self.mod for b in x)

class ImageUtils:
    """
    Clase que permite alterar imágenes
    """
    def __init__(self, image):
        """
        ImageUtils#ImageUtils(image : PIL.Image) : ImageUtils
        """
        self.image = image
    
    def negative(self):
        """
        ImageUtils#negative() : PIL.Image
        Devuelve el negativo de la imagen
        """
        # Invierte el valor de cada color de cada canal (R, G, B, posiblemente A) de la imagen
        neg = Image.eval(self.image, lambda c: 255 - c)
        
        # Si la imagen tiene un canal alpha de transparencia, recuperamos el original
        if self.image.mode in ('RGBA', 'LA'):
            neg.putalpha(self.image.split()[-1])
            
        return neg

##########################################################
# Práctica 1. Criptología                          #cripto
##########################################################
    
class Encoder:
    """
    Clase abstracta que representa un codificador de mensajes
    """
    def encode(self, message, *keys):
        """
        Encoder#encode(message : A, *keys : B) : A
        Método a implementar en las clases derivadas
        """
        raise NotImplementedError

    def decode(self, message, *keys):
        """
        Encoder#decode(message : A, *keys : B) : A
        Método a implementar en las clases derivadas
        """
        raise NotImplementedError

class Affine(Encoder):
    """
    Codificador/decodificador de imágenes mediante codificación afín
    """
    def encode(self, image, coeff, trans):
        """
        Affine#encode(image : PIL.Image, coeff : int, trans : (int, int, int)) : PIL.Image
        Codifica una imagen de acuerdo al coeficiente multiplicativo y
        el vector de traslación dados
        """
        # Separamos los canales de la imagen (requisito: imagen sin capa
        # de transparencia)
        canalr, canalg, canalb = image.split()

        # Realizamos la transformación afín en cada canal
        resr = Image.eval(canalr, lambda c: (coeff * c + trans[0]) % 256)
        resg = Image.eval(canalg, lambda c: (coeff * c + trans[1]) % 256)
        resb = Image.eval(canalb, lambda c: (coeff * c + trans[2]) % 256)

        # Volvemos a unir la imagen
        return Image.merge("RGB", (resr, resg, resb))

    def decode(self, image, coeff, trans):
        """
        Affine#decode(image : PIL.Image, coeff : int, trans : (int, int, int)) : PIL.Image
        Decodifica una imagen de acuerdo al coeficiente multiplicativo y
        el vector de traslación dados
        """
        # Separamos los canales de la imagen (requisito: imagen sin capa
        # de transparencia)
        canalr, canalg, canalb = image.split()

        # Calculamos el inverso multiplicativo del coeficiente en Z_256
        inv, _, _ = sy.gcdex(coeff, 256)

        # Deshacemos la transformación afín en cada canal
        resr = Image.eval(canalr, lambda c: (inv * (c - trans[0])) % 256)
        resg = Image.eval(canalg, lambda c: (inv * (c - trans[1])) % 256)
        resb = Image.eval(canalb, lambda c: (inv * (c - trans[2])) % 256)

        # Volvemos a unir la imagen
        return Image.merge("RGB", (resr, resg, resb))

class Vigenere(Encoder):
    """
    Codificador/decodificador Vigenère
    """
    def cycle(self, word):
        """
        Vigenere#cycle(word : string) : generator[char]
        Generador auxiliar que permite ciclar por una cadena de caracteres
        """
        n = len(word)
        i = 0
        while True:
            yield word[i % n]
            i += 1
    
    def encode(self, message, letters, keyword):
        """
        Vigenere#encode(message : string, letters : List[List[char]], keyword : string) : string
        Codifica mediante Vigenère para una matriz y una palabra clave dadas
        """
        # Almacenamos la primera columna para recorrerla más rápidamente
        first_column = [row[0] for row in letters]
        # Ciclamos por la palabra clave
        rows = self.cycle(keyword)
        
        return "".join([
            # Para cada letra del mensaje, encontramos en la matriz su codificación
            letters[
                # en la fila indicada por la letra correspondiente de la palabra clave
                first_column.index(rows.next())
            ][
                # y la columna indicada por la letra actual del mensaje
                letters[0].index(let)
            ]
            for let in message])

    def decode(self, message, letters, keyword):
        """
        Vigenere#decode(message : string, letters : List[List[char]], keyword : string) : string
        Decodifica mediante Vigenère para una matriz y una palabra clave dadas
        """
        # Almacenamos la primera columna para recorrerla más rápidamente
        first_column = [row[0] for row in letters]
        # Ciclamos por la palabra clave
        rows = self.cycle(keyword)

        return "".join([
            # Para cada letra del mensaje codificado, buscamos
            letters[
                # la primera fila
                0
            ][
                # la columna correspondiente a la posición ocupada por la letra
                # actual del mensaje en la fila indicada por la palabra clave
                letters[first_column.index(rows.next())].index(let)
            ]
            for let in message])

    def decode_keyword(self, message, letters, encoded):
        """
        Vigenere#decode_keyword(message : string, letters : List[List[char]], encoded : string) : string
        Toma un mensaje y su codificación por Vigenère y deduce la palabra
        clave que se usó
        """
        # Almacenamos la primera columna para recorrerla más rápidamente
        first_column = [row[0] for row in letters]
        
        return "".join([
            # Para cada letra del mensaje y su codificación, buscamos en la
            # primera columna
            first_column[
                # la fila correspondiente a la codificación en la columna
                # indicada por la letra del mensaje original
                letters[letters[0].index(message[i])].index(encoded[i])
            ]
            for i in xrange(len(message))])

class RSA(Encoder):
    """
    Codificador/decodificador RSA
    """
    def encode(self, message, mod, exponent):
        """
        RSA#encode(message : int, mod : int, exponent : int) : int
        Codifica mediante RSA
        """
        return pow(message, exponent, mod)

    def decode(self, message, p, q, exponent):
        """
        RSA#decode(message : int, p : int, q : int, exponent : int) : int
        Encuentra la phi de Euler de la clave privada (p, q) y la utiliza para
        descifrar el mensaje
        """
        n = p * q
        phi = (p - 1) * (q - 1)
        
        # Calculamos el exponente necesario al que elevaremos el mensaje
        # codificado para obtener el original. Este es el inverso multiplicativo
        # de `exponent` respecto de `phi`.
        u, _, _ = sy.gcdex(exponent, phi)

        return pow(message, long(u) % phi, n)
