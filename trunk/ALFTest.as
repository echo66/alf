﻿package {				import flash.display.MovieClip;	import flash.media.Sound;	import flash.text.TextField;	import flash.media.SoundChannel;	import flash.utils.ByteArray;	import flash.events.Event;	import flash.events.ProgressEvent;	import flash.events.SampleDataEvent;	import flash.net.*;	import flash.utils.Endian;	import flash.events.IOErrorEvent;	import flash.geom.ColorTransform;		import flash.geom.Rectangle;	import flash.net.FileReference;	import flash.events.*;	import flash.utils.*;	import flash.utils.Timer;	import flash.events.TimerEvent;	import fl.controls.Slider;	import fl.controls.CheckBox;	import fl.controls.ComboBox;	import fl.events.SliderEvent;	import flash.filters.BevelFilter;	import flash.geom.ColorTransform;	import ALF;		public class ALFTest extends MovieClip{						private var str:String;		//String for audio filename		private var myALF:ALF;		//ALF object		private var INITIALIZED:Boolean = false;				//Benchmarking		var frameNum:uint = 0;		public var time1:Number;		public var time2:Number;		public var total:Number;				//Features		public var inten:Number, cent:Number, band:Number, roll:Number, flux:Number;		var centArr:Array, bandArr:Array, intenArr:Array, rollArr:Array, fluxArr:Array;				//For Harmonics		public var freqs:Array, mags:Array;		public var numHarmonics:int = 1;				//Buttons		public var reverbStatus:Boolean = false;		public var reverbActive:String = "off";		public var harmonicsActive:Boolean = false;		public var audioPlaying:Boolean = false;				//Display		var roomRect:Rectangle = new Rectangle(0, 0, 200, 200);				//Room Stuff		private const roomXwidth:Number = 8; //8 meters wide		private const roomYwidth:Number = 8;		private const roomZheight:Number = 3;		private var srcX:Number = 0;		private var srcY:Number = 0;		private const srcZ:Number = 2;		//make the src height always 2 meters high		private var micX:Number = 0;		private var micY:Number = 0;		private const micZ:Number = 2;		//mic height should be approximately 2 meters, ~6feet		private var echoStrength:Number = 0.5;				public var songArray:Array;				//button look and feel		private var upBevel:BevelFilter = new BevelFilter(1, 45, 0xFFFFFF, 1, 0x000000, 1, 2, 2, 1, 1);		private var overBevel:BevelFilter = new BevelFilter(5, 45, 0xFFFFFF, 1, 0x000000, 1, 20, 20, 1, 3, 'inner', false);		private var downBevel:BevelFilter = new BevelFilter(-1, 45, 0xFFFFFF, 1, 0x000000, 1, 2, 2, 1, 1);				public function ALFTest(){			songArray = new Array();			songArray.push('http://music.ece.drexel.edu/~rmigneco/audio/SA1.wav');			songArray.push('http://music.ece.drexel.edu/~rmigneco/audio/IronAndWine.wav');			songArray.push('http://music.ece.drexel.edu/~rmigneco/audio/FurElise.wav');			songArray.push('http://music.ece.drexel.edu/~rmigneco/audio/lowE.wav');			//init the combo box			songComboBox.addItem({label:"--make a selection below--"});			songComboBox.addItem({label:"speech"});			songComboBox.addItem({label:"acoustic"});			songComboBox.addItem({label:"piano"});			songComboBox.addItem({label:"guitar string"});			songComboBox.addEventListener(Event.CHANGE, songChangedHandler);						//getHarmonics stuff			incHarmBtn.addEventListener(MouseEvent.CLICK, increaseHarmonics);			decHarmBtn.addEventListener(MouseEvent.CLICK, decreaseHarmonics);			numHarmonicsText.text = numHarmonics.toString();			incHarmBtn.addEventListener(MouseEvent.MOUSE_DOWN, buttonDown);			incHarmBtn.addEventListener(MouseEvent.MOUSE_UP, buttonUp);			decHarmBtn.addEventListener(MouseEvent.MOUSE_DOWN, buttonDown);			decHarmBtn.addEventListener(MouseEvent.MOUSE_UP, buttonUp);						//labeling			brightnessBtn.label = "";			fluxBtn.label = "";			bandwidthBtn.label = "";			intensityBtn.label = "";			rolloffBtn.label = "";						brightnessBtn.selected = false;			fluxBtn.selected = false;			bandwidthBtn.selected = false;			intensityBtn.selected = false;			rolloffBtn.selected = false;						// Storage Arrays//			centArr = new Array();			bandArr = new Array();			intenArr = new Array();			rollArr = new Array();			fluxArr = new Array();		}		public function onFrame(event:Event):void{			//check for selected features			if(brightnessBtn.selected) {								cent = myALF.getBrightness();				centArr.push(cent);				brightnessVal.text = (cent.toFixed(2)).toString() + ' Hz';			}			if(intensityBtn.selected) {				inten = myALF.getIntensity();				intenArr.push(inten);				intensityVal.text = (inten.toFixed(2)).toString();			}			if(rolloffBtn.selected) {				roll = myALF.getRolloff();				rollArr.push(roll);				rolloffVal.text = (roll.toFixed(2)).toString() + ' Hz';			}			if(fluxBtn.selected) {				flux = myALF.getFlux();				fluxArr.push(flux);				fluxVal.text = (flux.toFixed(2)).toString();			}			if(bandwidthBtn.selected) {				band = myALF.getBandwidth();				bandArr.push(band);				bandwidthVal.text = (band.toFixed(2)).toString() + ' Hz';			}						scaleOff();									//turn off any lit pitch buttons			if(harmonicsActive) { 						//harmonics handling code				myALF.getHarmonics(numHarmonics);		//return the desired number of harmonics				mags = myALF.getHarmonicAmplitudes();				freqs = myALF.getHarmonicFrequencies();				findPitch(freqs);						//code to light up the pitches			}						if(reverbActive == "on") {				//if reverb is on, get the src and mic positions 					micX = (theMic.x/theRoom.x * roomXwidth);				micY = (theMic.y/theRoom.y * roomYwidth);				srcX = (theSpeaker.x/theRoom.x * roomXwidth);				srcY = (theSpeaker.y/theRoom.y * roomYwidth);			}			//reverbDemo is undocumented, but allows full access to the DATF parameters			myALF.reverbDemo(reverbActive, 1, echoStrength, roomXwidth, roomYwidth, roomZheight,						 srcX, srcY, srcZ,						 micX, micY, micZ);			if(frameNum ==26){				myALF.getHarmonics(numHarmonics);		//return the desired number of harmonics				mags = myALF.getHarmonicAmplitudes();				freqs = myALF.getHarmonicFrequencies();				findPitch(freqs);						//code to light up the pitches							}						frameNum++;		}				public function audioLoaded(event:Event):void{						loadText.text = 'Loading Audio.... Complete!';			statusText.text = loadText.text;			initButtons();			playButton.addEventListener(MouseEvent.CLICK, playHandler);			playButton.alpha = 1;		}						public function initButtons():void{							loadText.text = 'Loading Audio.... Complete!';			playButton.buttonMode = true;			stopButton.buttonMode = true;			pauseButton.buttonMode = true;						//room reverb handling below...			reverbButton.addEventListener(MouseEvent.CLICK, reverbHandler);			theRoom.addChild(theSpeaker); 			//make the icons children of the room			theRoom.addChild(theMic);						//add click handling			theSpeaker.addEventListener(MouseEvent.MOUSE_DOWN, iconSelected);			theMic.addEventListener(MouseEvent.MOUSE_DOWN, iconSelected);			theSlider.addEventListener(SliderEvent.CHANGE, sliderValueChanged);						//initial positioins			theSpeaker.x = 50; theSpeaker.y = 50;			theMic.x = 100; theMic.y = 100;						//slider info			theSlider.maximum = 100;			theSlider.snapInterval = 10;			theSlider.tickInterval = 10;			//end initialization code						//audio control appearance			playButton.alpha = .33;			stopButton.alpha = .33;			pauseButton.alpha = .33;						//for harmoincs			noteButton.addEventListener(MouseEvent.MOUSE_DOWN, harmonicsToggled);			scaleOff();						//printing stuff, will be removed soon...			printFeaturesButton.alpha = .3;		}				public function songChangedHandler(evt:Event):void {			var index:int = evt.target.selectedIndex;			trace('You have selected: ' + songArray[index-1]);						statusText.text = "in song change";			centArr = [];			bandArr = []; 			intenArr = []; 			rollArr = []; 			fluxArr = []; 								// Initialize the ALF			if(index > 0){				playButton.removeEventListener(MouseEvent.CLICK, playHandler);				stopButton.removeEventListener(MouseEvent.CLICK	, stopHandler);				pauseButton.removeEventListener(MouseEvent.CLICK, pauseHandler);				playButton.alpha = .3; stopButton.alpha = .3; pauseButton.alpha = .3;				str = songArray[index - 1];				if(INITIALIZED){										myALF.loadNewSong(str);				}else{					myALF = new ALF(str, 10, true);					myALF.addEventListener(myALF.FILE_LOADED, audioLoaded);						myALF.addEventListener(myALF.PROG_EVENT, loadProgress);					myALF.addEventListener(myALF.NEW_FRAME, onFrame);					myALF.addEventListener(myALF.FILE_COMPLETE, audioFinished);					myALF.addEventListener(myALF.URL_ERROR, loadError);					INITIALIZED = true;				}			}		}				public function playHandler(event:Event):void{			if(!audioPlaying) {				audioPlaying = true;				//audio control handling below...				playButton.removeEventListener(MouseEvent.CLICK, playHandler);				stopButton.addEventListener(MouseEvent.CLICK, stopHandler);				pauseButton.addEventListener(MouseEvent.CLICK, pauseHandler);				songComboBox.enabled = false;				printFeaturesButton.removeEventListener(MouseEvent.CLICK, printFeaturesHandler);				printFeaturesButton.alpha = .3;				playButton.alpha = .30;				pauseButton.alpha = 1;				stopButton.alpha = 1;				myALF.startAudio();			}		}				public function stopHandler(event:Event):void{			if(audioPlaying) {				playButton.addEventListener(MouseEvent.CLICK, playHandler);				playButton.alpha = 1;				stopButton.alpha = .3;				pauseButton.alpha = .3;				songComboBox.enabled = true;				audioPlaying = false;				myALF.stopAudio();			}		}				public function pauseHandler(event:Event):void{			if(audioPlaying) {				playButton.addEventListener(MouseEvent.CLICK, playHandler);				playButton.alpha = 1;				pauseButton.alpha = .3;				stopButton.alpha =  .3;				audioPlaying = false;				myALF.pauseAudio();			}		}				public function audioFinished(event:Event):void{			audioPlaying = false;			playButton.addEventListener(MouseEvent.CLICK, playHandler);			printFeaturesButton.addEventListener(MouseEvent.CLICK, printFeaturesHandler);			printFeaturesButton.alpha = 1;			songComboBox.enabled = true;			playButton.alpha = 1;			stopButton.alpha = .3;			pauseButton.alpha = .3;			trace('audioFinished'); 			trace('---------------------------------------');			frameNum = 0;		}				public function sliderValueChanged(evt:SliderEvent):void {			echoStrength = evt.value/100;		}				public function harmonicsToggled(evt:Event):void {			if(harmonicsActive) {				harmonicsActive = false;				noteButton.filters = [overBevel, upBevel];			}			else {				harmonicsActive = true;				noteButton.filters = [overBevel, downBevel];			}		}		public function reverbHandler(evt:Event):void {			if(reverbStatus) {				//if its on, turn it off				reverbStatus = false;				reverbActive = "off";				reverbButton.filters = [overBevel, upBevel];			} else {				//of its off, turn it on				reverbStatus = true;				reverbActive = "on";				reverbButton.filters = [overBevel, downBevel];			}		}				public function loadProgress(event:Event):void {			var orgWidth:Number = 302.9;			loadText.text = 'Loading Audio.... ' + ((myALF.loadProgress*100).toFixed(2)).toString() + '%';			loaderBar.width = orgWidth * myALF.loadProgress;		}				public function loadError(evt:Event):void {			statusText.text = "Error loading specified URL. Check network connection. ";		}				//speaker/mic-icon handling functions		public function iconSelected(evt:Event):void {			evt.target.addEventListener(MouseEvent.MOUSE_UP, iconReleased);			evt.target.startDrag(false, roomRect);		}		public function iconReleased(evt:Event):void { evt.target.stopDrag(); }				public function buttonDown(evt:MouseEvent):void {			evt.target.alpha = .33;		}		public function buttonUp(evt:MouseEvent):void {			evt.target.alpha = 1;		}				public function increaseHarmonics(evt:MouseEvent):void {			if(numHarmonics < 10) {				incHarmBtn.addEventListener(MouseEvent.MOUSE_DOWN, buttonDown);				incHarmBtn.addEventListener(MouseEvent.MOUSE_UP, buttonUp);				incHarmBtn.alpha = 1; decHarmBtn.alpha = 1;				numHarmonics = numHarmonics + 1;				numHarmonicsText.text = numHarmonics.toString();				if(numHarmonics == 10) {					incHarmBtn.alpha = .33;					incHarmBtn.removeEventListener(MouseEvent.MOUSE_DOWN, buttonDown);					incHarmBtn.removeEventListener(MouseEvent.MOUSE_UP, buttonUp);				}			}		}				public function decreaseHarmonics(evt:MouseEvent):void {			if(numHarmonics > 1) {				decHarmBtn.addEventListener(MouseEvent.MOUSE_DOWN, buttonDown);				decHarmBtn.addEventListener(MouseEvent.MOUSE_UP, buttonUp);				incHarmBtn.alpha = 1; decHarmBtn.alpha = 1;				numHarmonics = numHarmonics - 1;				numHarmonicsText.text = numHarmonics.toString();				if(numHarmonics == 1) {					decHarmBtn.alpha = .33;					decHarmBtn.removeEventListener(MouseEvent.MOUSE_DOWN, buttonDown);					decHarmBtn.removeEventListener(MouseEvent.MOUSE_UP, buttonUp);				}			}		}				private function findPitch(freqArr:Array):void {			var freq:Number;			var numSteps:int;			var maxDist:Number;			var freqDist:Number;						for(var i:int = 0; i < freqArr.length; i++){				freq = freqArr[i];				maxDist = 40000;				//condense frequency to a 1 8va range				while(freq < 82.41 || freq > 164.82) {					if(freq < 82.41){ freq = freq * 2;}					else{ freq = freq/2;}				}								for(var step:int = 0; step < 13; step++) {					freqDist = Math.abs(freq - 82.41*Math.pow(2, step/12));					if(freqDist < maxDist) {						maxDist = freqDist;						numSteps = step;					}				}								//now turn on the appropriate button				switch(numSteps){					case 0:		eButton.alpha = 1;								break;					case 1:		fButton.alpha = 1;								break;					case 2:		fSharpButton.alpha = 1;								break;					case 3:		gButton.alpha = 1;								break;					case 4:		gSharpButton.alpha = 1;								break;					case 5:		aButton.alpha = 1;								break;					case 6:		aSharpButton.alpha = 1;								break;					case 7:		bButton.alpha = 1;								break;					case 8:		cButton.alpha = 1;								break;					case 9:		cSharpButton.alpha = 1;								break;					case 10:	dButton.alpha = 1;								break;					case 11:	dSharpButton.alpha = 1;								break;					case 12: 	eButton.alpha = 1;								break;					default: 	scaleOff();				}							}		}		private function scaleOff():void {			aButton.alpha = .30;			aSharpButton.alpha = .30;			bButton.alpha = .30;			cButton.alpha = .30;			cSharpButton.alpha = .30;			dButton.alpha = .30;			dSharpButton.alpha = .30;			eButton.alpha = .30;			fButton.alpha = .30;			fSharpButton.alpha = .30;			gButton.alpha = .30;			gSharpButton.alpha = .30;		}				public function printFeaturesHandler(evt:Event):void {			var file:FileReference = new FileReference();			var featureArray:Array = new Array();						if(brightnessBtn.selected) {							featureArray.push('*Centroid* ' + arr2str(centArr) + '\n');			}			if(fluxBtn.selected) {				featureArray.push('*Flux* ' + arr2str(fluxArr) + '\n');			}			if(bandwidthBtn.selected) {				featureArray.push('*Bandwidth* ' + arr2str(bandArr) + '\n');			}			if(intensityBtn.selected) {				featureArray.push('*Intensity* ' + arr2str(intenArr) + '\n');			}			if(rolloffBtn.selected) {				featureArray.push('*Rolloff* ' + arr2str(rollArr) + '\n');			}						file.save(featureArray, "features.txt");						centArr.splice(0); fluxArr.splice(0);			bandArr.splice(0); intenArr.splice(0); rollArr.splice(0); featureArray.splice(0);		}				public function arr2str(myArray:Array):String{			var i:uint;			var output:String = "";			for(i = 0; i < myArray.length; i++){ output = output+ " " + myArray[i].toString();}			return output;		}	}}