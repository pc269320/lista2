package com.example.dyrka2

import kotlin.math.pow

/**
 * Class Polynomial creates an object, which is a polynomial, based on given mutable list of
 * coefficients. It can calculate the degree of the polynomial, return it using the overridden
 * toString() method, calculate the value of the polynomial for a double. The class also allows
 * for Polynomial objects to be added, multiplied and subtracted. Returns error if no coefficients.
 */
class Polynomial(private var coefficients: MutableList<Double>){
    init{
        require(coefficients.isNotEmpty()) {"Input your coefficients"}
        while(coefficients.size>1 && coefficients.last()==0.0){ //usuwamy zera dla najwyższego x
            coefficients.removeAt(coefficients.size-1)
        }
    }

    /**
     * Function returns the degree of the polynomial.
     * Input:
     * -
     * Return:
     * polynomial degree (Int)
     */
    fun degree(): Int{
        return coefficients.size-1
    }

    /**
     * Function overrides the basic toString() method and allows for printing the polynomial in a
     * readable manner.
     * Input:
     * -
     * Return:
     * polynomial using math notation (String)
     */
    override fun toString(): String{
        var allZero = true
        for (coefficient in coefficients){
            if (coefficient != 0.0){
                allZero = false
                break
            }
        }
        if (allZero) return "P(x) = 0"

        val parts = mutableListOf<String>()
        for (i in coefficients.indices.reversed()){
            val coefficient = coefficients[i]
            if (coefficient != 0.0){
                val part = when (i){
                    0 -> "$coefficient"
                    1 -> "$coefficient*x"
                    else -> "$coefficient*x^$i"
                }
                parts.add(part)
            }
        }
        return "P(x) = " + parts.joinToString(" + ")
    }

    /**
     * Function calculates the value of a polynomial when it's invoked with a Double.
     * Input:
     * x (Double) - value for which the polynomial is calculated
     * Return:
     * result (Double) - value of the polynomial for the given input
     */
    operator fun invoke(x: Double): Double{
        var result= 0.0
        for (i in coefficients.indices){
            result+= coefficients[i] * x.pow((i).toDouble())
        }
        return result
    }

    /**
     * Function allows adding two Polynomial objects.
     * Input:
     * other (Polynomial) - another polynomial to add to the original polynomial
     * Return:
     * Polynomial(result) - product of Polynomial addition
     */
    operator fun plus(other: Polynomial): Polynomial{
        val maxSize = maxOf(this.coefficients.size, other.coefficients.size)
        val result = MutableList(maxSize){0.0}

        for (i in 0 until maxSize){
            val a = this.coefficients.getOrElse(i) {0.0}
            val b = other.coefficients.getOrElse(i) {0.0}
            result[i] = a + b
        }
        return Polynomial(result)
    }

    /**
     * Function allows subtracting two Polynomial objects.
     * Input:
     * other (Polynomial) - another polynomial to subtract to the original polynomial
     * Return:
     * Polynomial(result) - product of Polynomial subtraction
     */
    operator fun minus(other: Polynomial): Polynomial{
        val maxSize = maxOf(this.coefficients.size, other.coefficients.size)
        val result = MutableList(maxSize) {0.0}

        for (i in 0 until maxSize){
            val a = this.coefficients.getOrElse(i) {0.0}
            val b = other.coefficients.getOrElse(i) {0.0}
            result[i] = a - b
        }
        return Polynomial(result)
    }

    /**
     * Function allows multiplication of two Polynomial objects.
     * Input:
     * other (Polynomial) - another polynomial to multiply by the original polynomial
     * Return:
     * Polynomial(result) - product of Polynomial multiplication
     */
    operator fun times(other: Polynomial): Polynomial{
        val resultDegree = this.degree() + other.degree()
        val result = MutableList(resultDegree+1) {0.0}

        for (i in this.coefficients.indices){
            for (j in other.coefficients.indices){
                result[i+j] += this.coefficients[i] * other.coefficients[j]
            }
        }
        return Polynomial(result)
    }
}

fun main(){
    var p1 = Polynomial(mutableListOf(1.0, 2.0, 3.0, 4.0))
    println("\nPolynomial with coefficients (from highest degree to lowest) of {4, 3, 2, 1} is $p1")
    check(p1.toString()=="P(x) = 4.0*x^3 + 3.0*x^2 + 2.0*x + 1.0"){"Polynomials don't match!"}
    check(p1.degree()==3){"Polynomial degree is wrong!"}

    var p2 = Polynomial(mutableListOf(5.0, 10.0, 0.0, 0.0))
    println("\nPolynomial with coefficients (from highest degree to lowest) of {0, 0, 10, 5} is $p2")
    check(p2.toString()=="P(x) = 10.0*x + 5.0"){"Polynomials don't match!"}
    check(p2.degree()==1){"Polynomial degree is wrong!"}

    var p3 = Polynomial(mutableListOf(0.0))
    println("\nPolynomial with coefficient of zero is $p3")
    check(p3.toString()=="P(x) = 0"){"Polynomials don't match!"}
    check(p3.degree()==0){"Polynomial degree is wrong!"}

    var p4 = p1 + p2
    check(p4.toString()=="P(x) = 4.0*x^3 + 3.0*x^2 + 12.0*x + 6.0"){"Polynomials don't match!"}
    check(p4.degree()==3){"Polynomial degree is wrong!"}
    var p5 = p1
    p5+=p2
    check(p5.toString()==p4.toString()){"Polynomials don't match!"}

    p4 = p1 - p2
    check(p4.toString()=="P(x) = 4.0*x^3 + 3.0*x^2 + -8.0*x + -4.0"){"Polynomials don't match!"}
    check(p4.degree()==3){"Polynomial degree is wrong!"}
    p5 = p1
    p5-=p2
    check(p5.toString()==p4.toString()){"Polynomials don't match!"}

    p4 = p1 * p2
    println(p4)
    check(p4.toString()=="P(x) = 40.0*x^4 + 50.0*x^3 + 35.0*x^2 + 20.0*x + 5.0"){"Polynomials don't match!"}
    check(p4.degree()==4){"Polynomial degree is wrong!"}
    p5 = p1
    p5*=p2
    check(p5.toString()==p4.toString()){"Polynomials don't match!"}

    try{
        val pNoCoefficients = Polynomial(mutableListOf())
    } catch (e: IllegalArgumentException){
        println("\n${e.message}")
    }
}