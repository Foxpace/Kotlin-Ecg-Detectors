package com.motionapps.kotlin_ecg_detectors.realtime

abstract class RealTimeDetector(val fs: Double) {
    open fun processSample(sample: Double): Int{
        return -1
    }
    open fun reset(){}
}