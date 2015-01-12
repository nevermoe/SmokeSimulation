This is a fluid simulation program. The smoke is simulated using Euler method (grid-based method), and rendered with volume ray casting.


1. Prerequisites:
	
		opengl
		glew
		glfw
		Eigen


2. How to run:
		
		1. make
		2. ./main

3. Controll:

		Mouse:
		
			1. Change angle of view with mouse left key, and zoom with middle key.
			2. Select the Light and drag to change light position.

		Keyboard:
		
			1. R to reset the scene.
			2. S to switch between rendering and none rendering mode.
			3. W to toggle slices outline on/off.
			4. ESC to quit.

4. Screenshots:

	![ScreenShot](https://raw.githubusercontent.com/nevermoe/SmokeSimulation/master/screenshots/screenshot1.png)
 	![ScreenShot](https://raw.githubusercontent.com/nevermoe/SmokeSimulation/master/screenshots/screenshot2.png)
 	![ScreenShot](https://raw.githubusercontent.com/nevermoe/SmokeSimulation/master/screenshots/no_rendering.png)
 	![ScreenShot](https://raw.githubusercontent.com/nevermoe/SmokeSimulation/master/screenshots/drawframe.png)

