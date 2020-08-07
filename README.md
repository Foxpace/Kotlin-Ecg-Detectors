
# ECG R wave detectors in Kotlin

**This library is Kotlin translation of [py-ecg-detectors](https://github.com/berndporr/py-ecg-detectors).**

**_This Kotlin library is not created with affiliation to authors of Python version created by Luis Howell and Bernd Porr._**

**_Citation DOI: [10.5281/zenodo.3353396](https://doi.org/10.5281/zenodo.3353396) and license: [GPL-3.0 License](https://github.com/berndporr/py-ecg-detectors) to their work. Great gratitude for their work._** 

# Changes:

* main concepts of the algorithms stayed the same, only minor changes due to the implementation of the Kotlin
* **implementation of the FIR, IIR filters and other utilities are handled by [jDSP library for digital signal processing in Java](https://github.com/psambit9791/jDSP)**
* results of the algorithms may vary, due to the rounding errors and implementation of the filters 
* **missing swtDetector** due to the lack of appropriate wavelet decomposition algorithm - maybe will be added later
* HRV variability and testing is not available

### Available detectors:
* Hamilton detector
* Christov detector
* Engelse and Zeelenberg detector
* Pan and Tompkins detector
* Two moving averages detector
* Matched detector 
* *Missing - Stationary wavelet transform detector*

*Citations can be found above every function in Detectors class*

# Installation:

### Gradle version:

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
	        implementation 'com.github.Creative-Motion-Apps:kotlin-ecg-detectors:0.0.1'
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
	    <groupId>com.github.Creative-Motion-Apps</groupId>
	    <artifactId>kotlin-ecg-detectors</artifactId>
	    <version>0.0.1</version>
	</dependency>
```

# Usage:

```kotlin
    val fs = 512.0 // sampling frequency
    val array: DoubleArray = doubleArrayOf() // ECG signal 
    
    val detectors = Detectors(fs)  // main object for all detectors
    val rWaves: Array<Int> = detectors.engzeeDetector(array) // use of the specific detector
```


# Third parties and thanks goes to:


* **[jDSP](https://github.com/psambit9791/jDSP)** **- library focused on signal processing, filtration and manipulation in Java. Functionality is similar to scipy-signal package for Python. Library is created by *Sambit Paul*.** **_DOI: [10.5281/zenodo.3951042](https://doi.org/10.5281/zenodo.3951042). Library is licensed under MIT license._**

# License:

Licensed under the same license as [py-ecg-detectors](https://github.com/berndporr/py-ecg-detectors) - GPL-3.0 License. For more information visit LICENSE file in depository.
