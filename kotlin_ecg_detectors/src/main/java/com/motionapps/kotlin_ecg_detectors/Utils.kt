package com.motionapps.kotlin_ecg_detectors

import com.github.psambit9791.jdsp.misc.UtilMethods
import com.github.psambit9791.jdsp.signal.Convolution
import java.util.*
import kotlin.collections.ArrayList
import kotlin.math.abs

/**
 * Implementation of py-ecg-detectors in kotlin.
 * Copyright (C) 2020 creativemotion.app
 * GPL GNU GENERAL PUBLIC LICENSE Version 3, 29 June 2007
 */

object Utils {

    fun searchForMaximumInRange(doubleArray: DoubleArray, from: Int, to: Int): Double {
        var value = Double.NEGATIVE_INFINITY

        var f = from
        var t = to
        if(f < 0) f = 0
        if(t >= doubleArray.size) t = doubleArray.size - 1

        for (i in f..t) {
            if (doubleArray[i] > value) {
                value = doubleArray[i]
            }
        }
        return value
    }

    fun searchForMaximumInRangeArray(doubleArray: ArrayList<Double>, from: Int, to: Int): Double {
        var value = Double.NEGATIVE_INFINITY

        var f = from
        var t = to
        if(f < 0) f = 0
        if(t >= doubleArray.size) t = doubleArray.size - 1

        for (i in f..t) {
            if (doubleArray[i] > value) {
                value = doubleArray[i]
            }
        }
        return value
    }

    fun positiveFirstDifference(array: DoubleArray): DoubleArray {
        val arrayDif = DoubleArray(array.size)
        for (i in 1 until array.size) {
            arrayDif[i - 1] = abs(array[i] - array[i - 1])
        }
        arrayDif[arrayDif.lastIndex] = arrayDif[arrayDif.lastIndex - 1]
        return arrayDif
    }

    fun movedDifference(array: DoubleArray, move: Int): DoubleArray {
        val arrayDif = DoubleArray(array.size)
        for(i in move until array.size){
            arrayDif[i] = array[i] - array[i-move]
        }
        return arrayDif
    }

    fun positiveMiddleDifference(array: DoubleArray): DoubleArray {
        val values: ArrayList<Double> = ArrayList()
        for (i in 1 until array.size - 1) {
            values.add(abs(array[i + 1] - array[i - 1]))
        }
        return values.toDoubleArray()
    }


    fun addingZeros(array: DoubleArray, begin: Int, end: Int) {
        var b = begin
        var e = end

        if(b < 0) b = 0
        if(e >= array.size) e = array.size - 1

        for (i in b..e) {
            array[i] = 0.0
        }
    }

    fun meanArray(linkedList: LinkedList<Double>): Double {
        return linkedList.sum() / linkedList.size
    }

    fun pow2Array(array: DoubleArray): DoubleArray {
        val values = DoubleArray(array.size)
        for(i in array.indices){
            values[i] = array[i] * array[i]
        }
        return values
    }

    fun absArray(array: DoubleArray): DoubleArray{
        val values = DoubleArray(array.size)
        for(i in array.indices){
            values[i] = abs(array[i])
        }
        return values

    }

    fun MWA(array: DoubleArray, windowSize: Int): DoubleArray{
        val mwa = DoubleArray(array.size)
        for( i in array.indices){
            val section: DoubleArray = if( i < windowSize){
                UtilMethods.splitByIndex(array, 0, i)
            }else{
                UtilMethods.splitByIndex(array, i-windowSize, i)
            }

            if(i != 0){
                mwa[i] = section.sum()/section.size
            }else{
                mwa[i] = array[i]
            }
        }

        return mwa

    }

    fun valuesOfIndexes(array: DoubleArray, indexes: IntArray): DoubleArray{
        val values = DoubleArray(indexes.size)
        var counter = 0
        for(i in indexes){
            values[counter++] = array[i]
        }
        return values
    }

    fun valuesOfIndexesList(array: ArrayList<Double>, indexes: IntArray): DoubleArray{
        val values = DoubleArray(indexes.size)
        var counter = 0
        for(i in indexes){
            values[counter++] = array[i]
        }
        return values
    }

    fun diffInt(arr: IntArray): DoubleArray {
        val sig = DoubleArray(arr.size - 1)
        for (i in sig.indices) {
            sig[i] = (arr[i + 1] - arr[i]).toDouble()
        }
        return sig
    }

    class Smoother(private val signal: DoubleArray, wSize: Int) {

        private val kernel: DoubleArray = DoubleArray(wSize)
        private lateinit var output: DoubleArray

        private fun setKernel() {
            val value = 1.0 / kernel.size
            Arrays.fill(kernel, value)
        }

        fun smoothSignal(): DoubleArray {
            val c = Convolution(signal, kernel)
            output = c.convolve("same")
            return output
        }

        init {
            require(wSize <= signal.size) { "Kernel cannot be greater than signal." }
            setKernel()
        }
    }

}