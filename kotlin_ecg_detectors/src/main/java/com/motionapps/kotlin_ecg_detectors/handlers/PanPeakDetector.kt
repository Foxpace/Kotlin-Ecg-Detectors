package com.motionapps.kotlin_ecg_detectors.handlers

import com.motionapps.kotlin_ecg_detectors.Utils
import java.util.*

class PanPeakDetector(private val fs: Double) : HelperInterface {

    private val minDistance = (0.25 * fs)
    private val signalPeaks: LinkedList<Int> = LinkedList()
    private val noisePeaks: LinkedList<Int> = LinkedList()

    private var SPKI = 0.0
    private var NPKI = 0.0

    private var thresholdI1 = 0.0
    private var thresholdI2 = 0.0

    private var RRMissed = 0

    private var index = 0
    private val indexes: LinkedList<Int> = LinkedList()

    private val missedPeaks: LinkedList<Int> = LinkedList()
    private val peaks: LinkedList<Int> = LinkedList()

    init {
        signalPeaks.add(0)
    }

    private val filteredSignal: ArrayList<Double> = ArrayList()

    override fun passValue(sample: Double): Double {
        filteredSignal.add(sample)
        val i = filteredSignal.size - 2

        if (i > 0 && i < filteredSignal.size - 1) {
            if (filteredSignal[i - 1] < filteredSignal[i] && filteredSignal[i + 1] < filteredSignal[i]) {
                val peak = i
                peaks.add(i)

                if (filteredSignal[peak] > thresholdI1 && (peak - signalPeaks.last) > 0.3 * fs) {
                    signalPeaks.add(peak)
                    indexes.add(index)
                    SPKI = 0.125 * filteredSignal[signalPeaks.last] + 0.875 * SPKI

                    if (RRMissed != 0) {
                        if (signalPeaks.last - signalPeaks[signalPeaks.size - 2] > RRMissed) {
                            val missedSectionPeaks = LinkedList(
                                peaks.subList(
                                    indexes[indexes.size - 2] + 1,
                                    indexes.last
                                )
                            )
                            val missedSectionPeaks2 = LinkedList<Int>()

                            for (missedPeak in missedSectionPeaks) {
                                if (missedPeak - signalPeaks[signalPeaks.size - 2] > minDistance &&
                                    signalPeaks.last - missedPeak > minDistance &&
                                    filteredSignal[missedPeak] > thresholdI2
                                ) {
                                    missedSectionPeaks2.add(missedPeak)
                                }
                            }

                            if (missedSectionPeaks2.size > 0) {
                                val o = Utils.valuesOfIndexesList(
                                    filteredSignal,
                                    missedSectionPeaks2.toIntArray()
                                )
                                o.maxOrNull()?.let {
                                    val missedPeak =
                                        missedSectionPeaks2[o.indexOfFirst { it == it }]
                                    missedPeaks.add(missedPeak)
                                    signalPeaks.add(signalPeaks.last)
                                    signalPeaks[signalPeaks.size - 2] = missedPeak
                                    callAtEnd()
                                    return missedPeak.toDouble()
                                }
                            }
                        }
                    }
                    callAtEnd()
                    return peak.toDouble()
                } else {
                    noisePeaks.add(peak)
                    NPKI = 0.125 * filteredSignal[noisePeaks.last] + 0.875 * NPKI
                }
                callAtEnd()
            }
        }
        return -1.0

    }

    override fun reset() {
        signalPeaks.clear()
        noisePeaks.clear()
        SPKI = 0.0
        NPKI = 0.0
        thresholdI1 = 0.0
        thresholdI2 = 0.0
        RRMissed = 0
        index = 0
        indexes.clear()
        missedPeaks.clear()
        peaks.clear()
    }

    private fun callAtEnd() {
        thresholdI1 = NPKI + 0.25 * (SPKI - NPKI)
        thresholdI2 = 0.5 * thresholdI1

        if (signalPeaks.size > 8) {
            val array = LinkedList(
                signalPeaks.subList(
                    signalPeaks.size - 9,
                    signalPeaks.size - 1
                )
            ).toIntArray()
            val RR = Utils.diffInt(array)
            val RRAve = (RR.sum() / RR.size).toInt()
            RRMissed = (1.66 * RRAve.toDouble()).toInt()

        }

        index += 1
    }
}

