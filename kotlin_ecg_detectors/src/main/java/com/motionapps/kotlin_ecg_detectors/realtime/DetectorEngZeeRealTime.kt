package com.motionapps.kotlin_ecg_detectors.realtime

import com.github.psambit9791.jdsp.misc.UtilMethods
import com.motionapps.kotlin_ecg_detectors.Utils
import com.motionapps.kotlin_ecg_detectors.handlers.Convolution
import com.motionapps.kotlin_ecg_detectors.handlers.Filters
import com.motionapps.kotlin_ecg_detectors.handlers.HelperInterface
import com.motionapps.kotlin_ecg_detectors.handlers.MoveDifference
import java.util.*
import kotlin.collections.ArrayList

/**
 * C. Zeelenberg, A single scan algorithm for QRS detection and
 * feature extraction, IEEE Comp. in Cardiology, vol. 6,
 * pp. 37-42, 1979 with modifications A. Lourenco, H. Silva,
 * P. Leite, R. Lourenco and A. Fred, “Real Time
 * Electrocardiogram Segmentation for Finger Based ECG
 * Biometrics”, BIOSIGNALS 2012, pp. 49-54, 2012.
 */

class DetectorEngZeeRealTime(fs: Double, threshold: Double = 0.0085): RealTimeDetector(fs) {


    private val bandStopFilter = Filters.initBandStopFilter(this.fs, 4, 48.0, 52.0)
    private val moveDifference: HelperInterface = MoveDifference(4)

    private val kernel = doubleArrayOf(1.0, 4.0, 6.0, 4.0, 1.0)
    private val conv: HelperInterface = Convolution(kernel)

    private val ms200 = (0.2 * fs).toInt()
    private val ms1200 = (1.2 * fs).toInt()
    private val ms160 = (0.16 * fs).toInt()
    private val negThreshold = (threshold * fs).toInt()

    private var M = 0.0
    private var newM5 = 0.0
    private val MM: LinkedList<Double> = LinkedList()
    private val Mslope = UtilMethods.reverse(UtilMethods.arange(0.6, 1.0, 0.4 / (ms1200 - ms200).toDouble()))

    private val QRS: LinkedList<Int> = LinkedList()
    private val rPeaks: LinkedList<Int> = LinkedList()

    private var waitCounter: Int = (fs*0.2).toInt()
    private var counter = 0

    private val thiList: LinkedList<Int> = LinkedList()
    private var thi = false
    private var thf = false

    private val unfilteredValues: ArrayList<Double> = ArrayList()
    private val filteredValues: ArrayList<Double> = ArrayList()

    private var skipped = 0
    private var ready = false

    override fun processSample(sample: Double): Int{
        if (!ready) {
            var filteredSample = bandStopFilter.filter(sample)
            filteredSample = moveDifference.passValue(filteredSample)

            if(filteredSample.isNaN()){
                skipped++
                return -1
            }

            filteredSample = conv.passValue(filteredSample)
            if(filteredSample.isNaN()){
                skipped++
                return -1
            }

            if(waitCounter-- != 0){
                skipped++
                return -1
            }

            ready = true
            return -1
        } else {

            var filteredSample = bandStopFilter.filter(sample)
            filteredSample = moveDifference.passValue(filteredSample)
            filteredSample = conv.passValue(filteredSample)

            unfilteredValues.add(sample)
            filteredValues.add(filteredSample)
            val i = filteredValues.size - 1

            if (filteredValues.size < 5 * this.fs) {
                M = 0.6 * Utils.searchForMaximumInRangeArray(filteredValues, 0, i)
                MM.add(M)
                if (MM.size > 5) {
                    MM.removeFirst()
                }
            } else if (QRS.isNotEmpty() && filteredValues.size < QRS.last + ms200) {
                newM5 = 0.6 * Utils.searchForMaximumInRangeArray(filteredValues, QRS.last, i)
                if (newM5 > 1.5 * MM.last) {
                    newM5 = 1.1 * MM.last
                }

            } else if (QRS.isNotEmpty() && i == QRS.last + ms200) {
                MM.add(newM5)
                if (MM.size > 5) {
                    MM.removeFirst()
                }
                M = Utils.meanArray(MM)

            } else if (QRS.isNotEmpty() && i > QRS.last + ms200 && i < QRS.last + ms1200) {
                M = Utils.meanArray(MM) * Mslope[i - (QRS.last + ms200)]
            } else if (QRS.isNotEmpty() && i > QRS.last + ms1200) {
                M = 0.6 * Utils.meanArray(MM)
            }

            if (QRS.isEmpty() && filteredValues[i] > M) {
                QRS.add(i)
                thiList.add(i)
                thi = true
            } else if (QRS.isNotEmpty() && i > QRS.last + ms200 && filteredValues[i] > M) {
                QRS.add(i)
                thiList.add(i)
                thi = true
            }

            if (thi && i < thiList.last + ms160) {
                if (filteredValues[i] < -M && filteredValues[i - 1] > -M) {
                    thf = true
                }

                if (thf && filteredValues[i] < -M) {
                    counter++

                } else if (filteredValues[i] > -M && thf) {
                    counter = 0
                    thi = false
                    thf = false
                }
            } else if (thi && i > thiList.last + ms160) {
                counter = 0
                thi = false
                thf = false
            }

            if (counter > negThreshold) {
                val unfilteredSection = unfilteredValues.subList(thiList.last - negThreshold, i)
                unfilteredSection.maxOrNull()?.let {
                    val rPeak = unfilteredSection.indexOf(it) + thiList.last - negThreshold
                    rPeaks.add(rPeak)
                    counter = 0
                    thi = false
                    thf = false
                    return rPeak + skipped
                }
            }

            return -1
        }
    }

    override fun reset() {
        moveDifference.reset()
        conv.reset()
        M = 0.0
        newM5 = 0.0
        MM.clear()

        QRS.clear()
        rPeaks.clear()

        waitCounter = 0
        counter = 0

        thiList.clear()
        thi = false
        thf = false

        unfilteredValues.clear()
        filteredValues.clear()

        skipped = 0
        ready = false
    }

}