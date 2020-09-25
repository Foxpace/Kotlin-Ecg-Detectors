package com.motionapps.kotlin_ecg_detectors.handlers

import uk.me.berndporr.iirj.Butterworth
import kotlin.math.abs

object Filters {

    fun initBandPassFilter(
        samplingFrequency: Double,
        order: Int,
        lowCutoff: Double,
        highCutoff: Double): Butterworth {
        val centreFreq: Double = (highCutoff + lowCutoff) / 2.0
        val width: Double = abs(highCutoff - lowCutoff)
        val bp = Butterworth()
        bp.bandPass(order, samplingFrequency, centreFreq, width)
        return bp
    }

    fun initBandStopFilter(
        samplingFrequency: Double,
        order: Int,
        lowCutoff: Double,
        highCutoff: Double
    ): Butterworth {
        val centreFreq = (highCutoff + lowCutoff) / 2.0
        val width = abs(highCutoff - lowCutoff)
        val bs = Butterworth()
        bs.bandStop(order, samplingFrequency, centreFreq, width)
        return bs
    }

}