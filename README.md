
# ECG R wave detectors in Kotlin

**This library is Kotlin translation of [py-ecg-detectors](https://github.com/berndporr/py-ecg-detectors).**

**_This Kotlin library is not created with affiliation to authors of Python version created by Luis Howell and Bernd Porr._** 

**_Citation DOI: [10.5281/zenodo.3353396](https://doi.org/10.5281/zenodo.3353396) and license: [GPL-3.0 License](https://github.com/berndporr/py-ecg-detectors) to their work.  Great gratitude for their work._** [![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

# Changes:

* main concepts of the algorithms stayed the same, only minor changes due to the implementation of the Kotlin
* **implementation of the FIR, IIR filters and other utilities are handled by [jDSP library for digital signal processing in Java](https://github.com/psambit9791/jDSP)**
* results of the algorithms may vary, due to the rounding errors and implementation of the filters 
* **missing swtDetector** due to the lack of appropriate wavelet decomposition algorithm - maybe will be added later
* HRV variability and testing is not available
* real-time detectors use the same procedures customised to work with real-time samples / in for cycle

### Available offline detectors:
* Hamilton detector
* Christov detector
* Engelse and Zeelenberg detector
* Pan and Tompkins detector
* Two moving averages detector
* Matched detector 
* *Missing - Stationary wavelet transform detector*

*Citations can be found above every function in Detectors class*

### Available real-time detectors:
* Engelse and Zeelenberg detector
* Pan and Tompkins detector 

### [For C++ implementation go here](https://github.com/Creative-Motion-Apps/Cpp_Ecg_Detectors) 

# Installation:

### Gradle version:

Use this tag -> [![](https://jitpack.io/v/Foxpace/Kotlin-Ecg-Detectors.svg)](https://jitpack.io/#Foxpace/Kotlin-Ecg-Detectors)


**Jitpack.io dependency to 'build.gradle'**

```groovy
	allprojects {
		repositories {
			...
			maven { url 'https://jitpack.io' }
		}
	}

```

```groovy
	dependencies {
	        implementation 'com.github.Foxpace:Kotlin-Ecg-Detectors:Tag'
	}

```

### Maven version:

```xml
	<repositories>
		<repository>
		    <id>jitpack.io</id>
		    <url>https://jitpack.io</url>
		</repository>
	</repositories>
```

```xml
	<dependency>
	    <groupId>com.github.Foxpace</groupId>
	    <artifactId>Kotlin-Ecg-Detectors</artifactId>
	    <version>Tag</version>
	</dependency>
```

# Offline detector usage:

```kotlin
    val fs = 512.0 // sampling frequency
    val array: DoubleArray = doubleArrayOf() // ECG signal 
    
    val detectors = Detectors(fs)  // main object for all detectors
    val rWaves: Array<Int> = detectors.engzeeDetector(array) // use of the specific detector
```

# Real-time detector usage:
```kotlin
    val fs = 512.0 // sampling frequency
    val detector = DetectorEngZeeRealTime(fs, ...) // create specific detector with inputs
    val array: DoubleArray = doubleArrayOf() // ECG signal - just for demonstration
    // to simulate real-time procedure, we will use for
    for(sample in array){
        val qrs = detector.processSample(sample)
        if(qrs != -1){
            // -1 - nothing
            // otherwise - it is index of the QRS complex in signal
            // does not need to be in order for some detectors (sorting required in some cases)
        }
    } 
    detector.reset() // if something went wrong - this will empty memory in detector
    // for new configuration create new detector object
```

# Third parties and thanks goes to:


* **[jDSP](https://github.com/psambit9791/jDSP)** **- library focused on signal processing, filtration and manipulation in Java. Functionality is similar to scipy-signal package for Python. Library is created by *Sambit Paul*.** **_DOI: [10.5281/zenodo.3951042](https://doi.org/10.5281/zenodo.3951042). Library is licensed under MIT license._** [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
* **[iirj](https://github.com/berndporr/iirj)** **- IIR filters written in Java, which is also used by jDSP. Library is created by *Bernd Porr*, who is also contributor to py-ecg-detectors**. **_Library is licensed under Apache-2.0 License_**. [![License](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)

# License:

Licensed under the same license as [py-ecg-detectors](https://github.com/berndporr/py-ecg-detectors) - GPL-3.0 License. For more information visit LICENSE file in depository. [![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
