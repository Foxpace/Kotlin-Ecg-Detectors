package com.motionapps.kotlin_ecg_detectors.handlers

interface HelperInterface {
    fun passValue(sample: Double): Double
    fun reset()
}