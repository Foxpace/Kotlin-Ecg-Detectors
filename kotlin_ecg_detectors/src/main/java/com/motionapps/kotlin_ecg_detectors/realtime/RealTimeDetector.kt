package com.motionapps.kotlin_ecg_detectors.realtime

abstract class RealTimeDetector(val samplingFrequency: Double) {
    open fun processSample(sample: Double): Int{
        return -1
    }
}