package com.motionapps.kotlin_ecg_detectors.handlers

import org.apache.commons.math3.util.MathArrays

class Convolution(private val kernel: DoubleArray): HelperInterface {

    private val values: ArrayList<Double> = ArrayList()

    override fun passValue(sample: Double): Double {
        values.add(sample)
        return if (values.size == kernel.size) {
            val result = MathArrays.convolve(values.toDoubleArray(), kernel)[kernel.size - 1]
            values.removeFirst()
            result
        } else {
            Double.NaN
        }
    }

    override fun reset() {
        values.clear()
    }
}