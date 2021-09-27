package com.motionapps.kotlin_ecg_detectors

import com.github.psambit9791.jdsp.filter.Butterworth
import com.github.psambit9791.jdsp.misc.UtilMethods
import com.github.psambit9791.jdsp.signal.Convolution
import com.motionapps.kotlin_ecg_detectors.Utils.MWA
import com.motionapps.kotlin_ecg_detectors.Utils.absArray
import com.motionapps.kotlin_ecg_detectors.Utils.addingZeros
import com.motionapps.kotlin_ecg_detectors.Utils.diffInt
import com.motionapps.kotlin_ecg_detectors.Utils.meanArray
import com.motionapps.kotlin_ecg_detectors.Utils.movedDifference
import com.motionapps.kotlin_ecg_detectors.Utils.positiveFirstDifference
import com.motionapps.kotlin_ecg_detectors.Utils.positiveMiddleDifference
import com.motionapps.kotlin_ecg_detectors.Utils.pow2Array
import com.motionapps.kotlin_ecg_detectors.Utils.searchForMaximumInRange
import com.motionapps.kotlin_ecg_detectors.Utils.valuesOfIndexes
import java.util.*
import kotlin.math.roundToInt

/* License from: https://github.com/berndporr/py-ecg-detectors/blob/master/LICENSE */
/**
* A collection of 7 ECG heartbeat detection algorithms implemented
* in Python. Developed in conjunction with a new ECG database:
* http://researchdata.gla.ac.uk/716/.
* Copyright (C) 2019 Luis Howell & Bernd Porr
* GPL GNU GENERAL PUBLIC LICENSE Version 3, 29 June 2007
*/


class Detectors(private val fs: Double) {

    /**
     * P.S. Hamilton,
     * Open Source ECG Analysis Software Documentation, E.P.Limited, 2002.
     */

    fun hamiltonDetector(unfilteredSignal: DoubleArray): Array<Int> {

        // 8 Hz LP and 16 Hz HP filter
        val butterFilter = Butterworth(unfilteredSignal, fs)
        var filteredSignal: DoubleArray = butterFilter.bandPassFilter(1, 8.0, 16.0)
        filteredSignal = positiveFirstDifference(filteredSignal)

        // smoothing signal with 80 ms latency
        val sizeOfSmoothing = (0.08 * fs).toInt()
        Utils.Smoother(filteredSignal, sizeOfSmoothing).also {
            filteredSignal = it.smoothSignal()
        }

        addingZeros(filteredSignal, begin = 0, end = sizeOfSmoothing * 2)

        val nPks: LinkedList<Double> = LinkedList()
        var nPksAve = 0.0 // average of noise peaks - peaks which are not classified as QRS

        val sPks: LinkedList<Double> = LinkedList()
        var sPksAve = 0.0 // average of recent QRS peaks

        val QRS: LinkedList<Int> = LinkedList()
        QRS.add(0)
        val RR: LinkedList<Int> = LinkedList()
        var RRAve = 0.0

        var threshold = 0.0

        val idx: LinkedList<Int> = LinkedList()
        val peaks: LinkedList<Int> = LinkedList()

        for (peak in 1 until filteredSignal.size) {
            // peak detection
            if (filteredSignal[peak - 1] < filteredSignal[peak] && filteredSignal[peak] > filteredSignal[peak + 1]) {
                peaks.add(peak)

                //  signal needs to be above threshold, but peak under 360 ms from previous detection is taken as T wave
                if (filteredSignal[peak] > threshold && (peak - QRS.last) > 0.3 * fs) {
                    QRS.add(peak)
                    idx.add(peak)
                    sPks.add(filteredSignal[peak])

                    if (nPks.size > 8) {
                        sPks.removeFirst()
                    }

                    sPksAve = meanArray(sPks)

                    if (RRAve != 0.0) {
                        // if in interval 1.5 R-R, not R detection, search for peak with halved threshold
                        if (QRS.last - QRS[QRS.size - 2] > 1.5 * RRAve) {
                            try {
                                val missedPeaks: MutableList<Int> = peaks.subList(idx[idx.size - 2] + 1, idx.last)
                                for (missedPeak in missedPeaks) {
                                    if (missedPeak - peaks[idx[idx.size - 2]] > (0.36 * fs).toInt() && filteredSignal[missedPeak] > 0.5 * threshold) {
                                        QRS.add(missedPeak)
                                        QRS.sort()
                                        break
                                    }
                                }
                            }catch (e: IndexOutOfBoundsException){}
                        }
                    }

                    if (QRS.size > 2) {
                        RR.add(QRS.last - QRS[QRS.size - 2])
                        if (RR.size > 8) {
                            RR.removeFirst()
                        }
                        RRAve = RR.sum().toDouble() / RR.size.toDouble()
                    }

                } else {
                    nPks.add(filteredSignal[peak])
                    if (nPks.size > 8) {
                        nPks.removeFirst()
                    }
                    nPksAve = meanArray(nPks)
                }

                // threshold = average_noise_peak + 0.45 (average_QRS_peak - average_noise_peak)
                threshold = nPksAve + 0.45 * (sPksAve - nPksAve)

            }
        }

        QRS.removeFirst()

        return QRS.toTypedArray()

    }

    /**
     * Ivaylo I. Christov,
     * Real time electrocardiogram QRS detection using combined
     * adaptive threshold, BioMedical Engineering OnLine 2004,
     * vol. 3:28, 2004.
     */

    fun christovDetector(unfilteredSignal: DoubleArray): Array<Int> {

        var totalTaps = 0

        var sizeOfSmoothing = (0.02 * fs).toInt()
        var smoother: Utils.Smoother = Utils.Smoother(unfilteredSignal, sizeOfSmoothing)
        totalTaps += sizeOfSmoothing

        val MA1: DoubleArray = smoother.smoothSignal()

        sizeOfSmoothing = (0.028 * fs).toInt()
        smoother = Utils.Smoother(MA1, sizeOfSmoothing)
        totalTaps += sizeOfSmoothing

        val MA2: DoubleArray = smoother.smoothSignal()

        val Y = positiveMiddleDifference(MA2)

        sizeOfSmoothing = (0.04 * fs).toInt()
        smoother = Utils.Smoother(Y, sizeOfSmoothing)
        totalTaps += sizeOfSmoothing

        val MA3: DoubleArray = smoother.smoothSignal()

        addingZeros(MA3, 0, totalTaps)

        val ms50 = (0.05 * fs).toInt()
        val ms200 = (0.2 * fs).toInt()
        val ms1200 = (1.2 * fs).toInt()
        val ms350 = (0.35 * fs).toInt()

        var M = 0.0
        var newM5 = 0.0
        val MM: LinkedList<Double> = LinkedList()
        val Mslope = UtilMethods.reverse(UtilMethods.arange(0.6, 1.0, 0.4 / (ms1200 - ms200).toDouble()))

        var F = 0.0

        var R = 0.0
        val RR: LinkedList<Int> = LinkedList()
        var RM = 0

        val QRS: LinkedList<Int> = LinkedList()

        for (i in MA3.indices) {

            if (i < 5 * fs) {

                M = 0.6 * searchForMaximumInRange(MA3, 0, i + 1)
                MM.add(M)
                if (MM.size > 5) {
                    MM.removeFirst()
                }

            } else if (QRS.isNotEmpty() && i < (QRS.last + ms200)) {
                newM5 = 0.6 * searchForMaximumInRange(MA3, QRS.last, i)
                if (newM5 > 1.5 * MM.last) {
                    newM5 = 1.1 * MM.last
                }
            } else if (QRS.isNotEmpty() && i == (QRS.last + ms200)) {
                if (newM5 == 0.0) {
                    newM5 = MM.last
                }
                MM.add(newM5)
                if (MM.size > 5) {
                    MM.removeFirst()
                }

                M = meanArray(MM)

            } else if (QRS.isNotEmpty() && i > (QRS.last + ms200) && i < (QRS.last + ms1200)) {
                M = meanArray(MM) * Mslope[i - (QRS.last + ms200)]
            } else if (QRS.isNotEmpty() && i > (QRS.last + ms1200)) {
                M = 0.6 * meanArray(MM)
            }

            if (i > ms350) {
                val fSection = MA3.copyOfRange(i - ms350, i)
                val maxLatest = searchForMaximumInRange(fSection, fSection.size - ms50, fSection.size - 1)
                val maxEarliest = searchForMaximumInRange(fSection, 0, ms50)
                F += (maxLatest - maxEarliest) / 150.0
            }

            if (QRS.isNotEmpty() && i < (QRS.last + ((2.0 / 3.0) * RM).toInt())) {
                R = 0.0
            } else if (QRS.isNotEmpty() && i > (QRS.last + ((2.0 / 3.0) * RM).toInt()) && i < (QRS.last + RM)) {
                R = 0 + (M - meanArray(MM)) / 1.4
            }

            val MFR = M + F + R

            if (QRS.isEmpty() && (MA3[i] > MFR)) {
                QRS.add(i)

            } else if (QRS.isNotEmpty() && i > (QRS.last + ms200) && (MA3[i] > MFR)) {

                QRS.add(i)
                if (QRS.size > 2) {
                    RR.add(QRS.last - QRS[QRS.size - 2])
                    if (RR.size > 5) {
                        RR.removeFirst()
                    }
                    RM = (RR.sum().toDouble() / RR.size.toDouble()).roundToInt()
                 }
            }
        }

        QRS.removeFirst()

        return QRS.toTypedArray()
    }

    /**
     * C. Zeelenberg, A single scan algorithm for QRS detection and
     * feature extraction, IEEE Comp. in Cardiology, vol. 6,
     * pp. 37-42, 1979 with modifications A. Lourenco, H. Silva,
     * P. Leite, R. Lourenco and A. Fred, “Real Time
     * Electrocardiogram Segmentation for Finger Based ECG
     * Biometrics”, BIOSIGNALS 2012, pp. 49-54, 2012.
     */

    fun engzeeDetector(unfilteredSignal: DoubleArray, threshold: Double = 0.0085): Array<Int> {


        // band-stop 48 - 52 Hz
        val butter = Butterworth(unfilteredSignal, fs)
        var filteredSignal = butter.bandStopFilter(4, 48.0, 52.0)

        filteredSignal = movedDifference(filteredSignal, 4)

        val ci = doubleArrayOf(1.0, 4.0, 6.0, 4.0, 1.0)
        Convolution(filteredSignal, ci).also {
            filteredSignal = it.convolve("same")
        }

        addingZeros(filteredSignal, 0, (0.2 * fs).toInt())

        val ms200 = (0.2 * fs).toInt()
        val ms1200 = (1.2 * fs).toInt()
        val ms160 = (0.16 * fs).toInt()
        val negThreshold = (threshold * fs).toInt()

        var M = 0.0
        var newM5 = 0.0
        val MM: LinkedList<Double> = LinkedList()
        val Mslope =
            UtilMethods.reverse(UtilMethods.arange(0.6, 1.0, 0.4 / (ms1200 - ms200).toDouble()))

        val QRS: LinkedList<Int> = LinkedList()
        val rPeaks: LinkedList<Int> = LinkedList()

        var counter = 0

        val thiList: LinkedList<Int> = LinkedList()
        var thi = false
        var thf = false

        for (i in filteredSignal.indices) {

            if (i < 5 * fs) {
                M = 0.6 * searchForMaximumInRange(filteredSignal, 0, i)
                MM.add(M)
                if (MM.size > 5) {
                    MM.removeFirst()
                }
            } else if (QRS.isNotEmpty() && i < QRS.last + ms200) {
                newM5 = 0.6 * searchForMaximumInRange(filteredSignal, QRS.last, i)
                if (newM5 > 1.5 * MM.last) {
                    newM5 = 1.1 * MM.last
                }

            } else if (QRS.isNotEmpty() && i == QRS.last + ms200) {
                MM.add(newM5)
                if (MM.size > 5) {
                    MM.removeFirst()
                }
                M = meanArray(MM)

            } else if (QRS.isNotEmpty() && i > QRS.last + ms200 && i < QRS.last + ms1200) {
                M = meanArray(MM) * Mslope[i - (QRS.last + ms200)]
            } else if (QRS.isNotEmpty() && i > QRS.last + ms1200) {
                M = 0.6 * meanArray(MM)
            }

            if (QRS.isEmpty() && filteredSignal[i] > M) {
                QRS.add(i)
                thiList.add(i)
                thi = true
            } else if (QRS.isNotEmpty() && i > QRS.last + ms200 && filteredSignal[i] > M) {
                QRS.add(i)
                thiList.add(i)
                thi = true
            }

            if (thi && i < thiList.last + ms160) {
                if (filteredSignal[i] < -M && filteredSignal[i-1] > -M) {
                    thf = true
                }

                if (thf && filteredSignal[i] < -M) {
                    counter++

                } else if (filteredSignal[i] > -M && thf) {
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
                val unfilteredSection = unfilteredSignal.copyOfRange(thiList.last - negThreshold, i)
                unfilteredSection.maxOrNull()?.let {
                    rPeaks.add(unfilteredSection.indexOfFirst { it == it } + thiList.last - negThreshold + (0.03 * fs).toInt())
                }

                counter = 0
                thi = false
                thf = false
            }
        }

        return rPeaks.toTypedArray()
    }

    /**
     * Jiapu Pan and Willis J. Tompkins.
     * A Real-Time QRS Detection Algorithm.
     * In: IEEE Transactions on Biomedical Engineering
     * BME-32.3 (1985), pp. 230–236.
     */

    fun panTompkinsDetector(unfilteredSignal: DoubleArray): Array<Int> {

        val butter = Butterworth(unfilteredSignal, fs)
        var filteredSignal = butter.bandPassFilter(1, 5.0, 15.0)
        filteredSignal = UtilMethods.diff(filteredSignal)

        filteredSignal = pow2Array(filteredSignal)

        val N = (0.12 * fs).toInt()
        val mwa = MWA(filteredSignal, N)
        addingZeros(mwa, 0, (0.2 * fs).toInt())

        return panPeakDetect(mwa, fs)


    }

    private fun panPeakDetect(detection: DoubleArray, fs: Double): Array<Int> {
        val minDistance = (0.25 * fs)

        val signalPeaks: LinkedList<Int> = LinkedList()
        signalPeaks.add(0)

        val noisePeaks: LinkedList<Int> = LinkedList()

        var SPKI = 0.0
        var NPKI = 0.0

        var thresholdI1 = 0.0
        var thresholdI2 = 0.0

        var RRMissed = 0
        var index = 0
        val indexes: LinkedList<Int> = LinkedList()

        val missedPeaks: LinkedList<Int> = LinkedList()
        val peaks: LinkedList<Int> = LinkedList()

        for (i in detection.indices) {

            if (i > 0 && i < detection.size - 1) {
                if (detection[i-1] < detection[i] && detection[i + 1] < detection[i]) {
                    val peak = i
                    peaks.add(i)

                    if (detection[peak] > thresholdI1 && (peak - signalPeaks.last) > 0.3 * fs) {
                        signalPeaks.add(peak)
                        indexes.add(index)
                        SPKI = 0.125 * detection[signalPeaks.last] + 0.875 * SPKI

                        if (RRMissed != 0) {
                            if (signalPeaks.last - signalPeaks[signalPeaks.size - 2] > RRMissed) {
                                val missedSectionPeaks = LinkedList(
                                    peaks.subList(
                                        indexes[indexes.size - 2]+1,
                                        indexes.last
                                    )
                                )
                                val missedSectionPeaks2 = LinkedList<Int>()

                                for (missedPeak in missedSectionPeaks) {
                                    if (missedPeak - signalPeaks[signalPeaks.size - 2] > minDistance &&
                                        signalPeaks.last - missedPeak > minDistance &&
                                        detection[missedPeak] > thresholdI2
                                    ) {
                                        missedSectionPeaks2.add(missedPeak)
                                    }
                                }

                                if (missedSectionPeaks2.size > 0) {
                                    val o = valuesOfIndexes(detection, missedSectionPeaks2.toIntArray())
                                    o.maxOrNull()?.let {
                                        val missedPeak = missedSectionPeaks2[o.indexOfFirst { it == it }]
                                        missedPeaks.add(missedPeak)
                                        signalPeaks.add(signalPeaks.last)
                                        signalPeaks[signalPeaks.size - 2] = missedPeak
                                    }
                                }
                            }
                        }
                    } else {
                        noisePeaks.add(peak)
                        NPKI = 0.125 * detection[noisePeaks.last] + 0.875 * NPKI

                    }

                    thresholdI1 = NPKI + 0.25 * (SPKI - NPKI)
                    thresholdI2 = 0.5 * thresholdI1

                    if (signalPeaks.size > 8) {
                        val array = LinkedList(
                            signalPeaks.subList(
                                signalPeaks.size - 9,
                                signalPeaks.size - 1
                            )
                        ).toIntArray()
                        val RR = diffInt(array)
                        val RRAve = (RR.sum() / RR.size).toInt()
                        RRMissed = (1.66 * RRAve.toDouble()).toInt()

                    }

                    index += 1

                }
            }

        }
        signalPeaks.removeFirst()
        return signalPeaks.toTypedArray()
    }

    /**
     * Elgendi, Mohamed & Jonkman,
     * Mirjam & De Boer, Friso. (2010).
     * Frequency Bands Effects on QRS Detection.
     * The 3rd International Conference on Bio-inspired Systems
     * and Signal Processing (BIOSIGNALS2010). 428-431.
     * */
    fun twoAverageDetector(unfilteredSignal: DoubleArray): Array<Int>{

        val butter = Butterworth(unfilteredSignal, fs)
        val filteredSignal = butter.bandPassFilter(2, 8.0, 20.0)

        val window1 = (0.12*fs).toInt()
        val mwaQrs = MWA(absArray(filteredSignal), window1)

        val window2 = (0.6*fs).toInt()
        val mwaBeat = MWA(absArray(filteredSignal), window2)

        val blocks = DoubleArray(unfilteredSignal.size)
        val blockHeight: Double = unfilteredSignal.maxOrNull()!!

        for (i in mwaQrs.indices){
            if(mwaQrs[i] > mwaBeat[i]){
                blocks[i] = blockHeight
            }else{
                blocks[i] = 0.0
            }
        }

        val QRS = LinkedList<Int>()

        var start = 0

        for(i in 1 until blocks.size){
            if(blocks[i-1] == 0.0 && blocks[i] == 0.0){
                start = i
            }else if(blocks[i-1] == blockHeight && blocks[i] == 0.0){
                val end = i-1
                if( end - start > (0.08*fs).toInt()){
                    val partOfSignal  = filteredSignal.copyOfRange(start, end+1)
                    partOfSignal.maxOrNull()?.let {
                        val detection = partOfSignal.indexOfFirst { it == it } +start
                        if(QRS.isNotEmpty()){
                            if(detection - QRS.last > (0.3*fs)){
                                QRS.add(detection)
                            }
                        }else{
                            QRS.add(detection)
                        }
                        true
                    }
                }
            }
        }

        return QRS.toTypedArray()

    }
    /**
     * Stationary Wavelet Transform
     * based on Vignesh Kalidas and Lakshman Tamil.
     * Real-time QRS detector using Stationary Wavelet Transform
     * for Automated ECG Analysis.
     * In: 2017 IEEE 17th International Conference on
     * Bioinformatics and Bioengineering (BIBE).
     * Uses the Pan and Tompkins threshold.
     * */

//    fun swtDetector(unfilteredSignal: DoubleArray): DoubleArray? {
//        val swtLevel = 3.0
//        var padding = (2.0.pow(ceil(log2(unfilteredSignal.size.toDouble()))) - unfilteredSignal.size.toDouble()).toInt()
//
//        val paddedSignal = if(padding > 0){
//            val cons = DoubleArray(padding)
//            Arrays.fill(cons, unfilteredSignal.last())
//            UtilMethods.concatenateArray(unfilteredSignal, cons)
//        }else{
//            unfilteredSignal
//        }
//        val transform = Transform(FastWaveletTransform(Daubechies3()))
//        var wfOutput = transform.forward(paddedSignal, swtLevel.toInt())
//
//        wfOutput = pow2Array(wfOutput)
//
//        val butter = Butterworth(wfOutput, fs)
//        wfOutput = butter.band_pass_filter(3, 0.01, 10.0)
//
////        return panPeakDetect(wfOutput, fs)
//
//        return wfOutput
//    }

    /**
     * Template of ECG wave is used for convolution with ECG signal. Peaks then represent QRS
     * presence, which are detected with Pan and Tompkins threshold method.
     * */

    fun matchedFilterDetector(unfilteredSignal: DoubleArray, template: DoubleArray): Array<Int>{

        val butter = Butterworth(unfilteredSignal, fs)
        val filteredSignal = butter.bandPassFilter(4, 0.1, 48.0)
        val reversedTemplate = UtilMethods.reverse(template)
        Convolution(filteredSignal, reversedTemplate).also {
            var detection = it.convolve("same")
            detection = pow2Array(detection)
            addingZeros(detection, 0, reversedTemplate.size)

            return panPeakDetect(detection, fs)
        }

    }

}