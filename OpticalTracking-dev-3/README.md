# Optical Tracking using Polaris Vega

## Usage
Take a look at `defaults.m`. These likely do not need to be changed, however you should look into the `Settings` and `Run flags` sections.

Add `lib/` to path, then run `main.m`.
A window will popup to the root folder where all specimen folders are.
Every specimen must include a folder `Digitisation` or `Calibration` with landmark that are reasonably identifiable: e.g., `tibia-medial`, `Tibia lateral`, `femur_proximal`, `fem_med`, `patella sup`, `TD`.

## My loading conditions are named inconsistently. What do I do?
in `lib/configure/clean_specimen_condition.m`, use `regexprep()` and `replace()` to correct names. Some useful examples:
```matlab
O = regexprep(O, 'ante[rior]*', 'Anterior', 'ignorecase'); % Correct various misspelings of the word "anterior", e.g. anteriro, anteioror, etc
O = regexprep(O, 'external .*', 'External', 'ignorecase'); % converts 'external' followed by any character into External. e.g., External Rot, external rotation, EXTERNAL Rota => External
```

## Overview
To understand what this code does, take a look at this image. We're assuming the z-axis is 0 to make life easier.

![rigid-bodies-demonstration](https://github.com/user-attachments/assets/995f4cdf-66b9-4d40-a1e5-de14d93a6faa)


The blue and red are objects in space with orientations based on their x and y axes.
The optical tracking camera gives us the black x,y axes (the _global_ frame of reference), but we are interested in the motion of the blue object relative to the red.

The origin of the ![blue](https://placehold.co/15x15/0000ff/0000ff.png) blue object could be described as `(2,1,0)` in the ![black](https://placehold.co/15x15/000/000.png)global reference, or `(1,-1, 0)` relative to the ![red](https://placehold.co/15x15/ff0000/ff0000.png) red object.

![rigid-bodies-demonstration-arrows](https://github.com/user-attachments/assets/a9bd3416-8e03-45f2-9baf-b62679a21f6f)


In our notation, these would be:
```math
_gr_b = \begin{bmatrix} 2 \\\ 1 \\\ 0 \end{bmatrix} \quad \quad _rr_b = \begin{bmatrix} 1 \\\ -1 \\\ 0 \end{bmatrix}
```

### What about rotation?
Let's consider ![red](https://placehold.co/15x15/ff0000/ff0000.png) red object first to show what no rotation looks like.
The transform that represents rotations and translations is a 4x4 matrix. 
```math
_gT_r = \begin{bmatrix} R & t \\\ 0 & 1\end{bmatrix} = 
\begin{bmatrix}
R_{xx} & R_{xy} & R_{xz} & t_x \\\
R_{yx} & R_{yy} & R_{yz} & t_y \\\
R_{zx} & R_{zy} & R_{zz} & t_z \\\
0 & 0 & 0 & 1
\end{bmatrix} 
```
`t`, the translations, are x,y,z translations from the global frame of reference to the object. From the image we can see that is `x: 1`, `y: 2`, `z: 0`
```math
_gT_r = 
\begin{bmatrix}
R_{xx} & R_{xy} & R_{xz} & 1 \\\
R_{yx} & R_{yy} & R_{yz} & 2 \\\
R_{zx} & R_{zy} & R_{zz} & 0 \\\
0 & 0 & 0 & 1
\end{bmatrix} 
```

each column of `R`, the rotations, shows how the unit vectors are orientated relative to the global frame of reference. We can see that the unit vectors for both the global and the red object are orientated the same way.
```math
\quad _g\hat{i}_g  = \begin{bmatrix} 1 \\\ 0 \\\ 0 \end{bmatrix} \quad\quad \text{and} \quad\quad _g\hat{i}_r = \begin{bmatrix} 1 \\\ 0 \\\ 0 \end{bmatrix} \quad\quad\quad\quad\\\
\quad _g\hat{j}_g  = \begin{bmatrix} 0 \\\ 1 \\\ 0 \end{bmatrix} \quad\quad \text{and} \quad\quad _g\hat{j}_r = \begin{bmatrix} 0 \\\ 1 \\\ 0 \end{bmatrix} \quad\text{...etc}
```
Which means we end up with a ![red](https://placehold.co/15x15/ff0000/ff0000.png) red transformation matrix that looks something like:
```math
_gT_r = 
\begin{bmatrix}
1 & 0 & 0 & 1 \\\
0 & 1 & 0 & 2 \\\
0 & 0 & 1 & 0 \\\
0 & 0 & 0 & 1
\end{bmatrix} 
```

If we apply the same logic to the blue object and assume, say, 60 degrees, then the `i` unit vector (x axis) has rotated like so:

![unit-vectors](https://github.com/user-attachments/assets/9ef4b804-decf-4440-99e8-ecb51a2953d5)


```math
 _g\hat{i}_b = \begin{bmatrix}
\frac{1}{2} \\\ \frac{\sqrt{3}}{2} \\\ 0
\end{bmatrix}
```
Doing the same to `j` and `k` and including the translation vector to create the transform that describes, in the ![black](https://placehold.co/15x15/000/000.png) global frame of reference, the ![blue](https://placehold.co/15x15/0000ff/0000ff.png) blue object:
```math
_gT_b = 
\begin{bmatrix}
\frac{1}{2} & -\frac{\sqrt{3}}{2}  & 0 & 2 \\\
\frac{\sqrt{3}}{2} & \frac{1}{2} & 0 & 1 \\\
0 & 0 & 1 & 0 \\\
0 & 0 & 0 & 1
\end{bmatrix} 
```
or, in the ![red](https://placehold.co/15x15/ff0000/ff0000.png) red frame of reference, the ![blue](https://placehold.co/15x15/0000ff/0000ff.png) blue object has this transform:
```math
_rT_b = 
\begin{bmatrix}
\frac{1}{2} & -\frac{\sqrt{3}}{2}  & 0 & 1 \\\
\frac{\sqrt{3}}{2} & \frac{1}{2} & 0 & -1 \\\
0 & 0 & 1 & 0 \\\
0 & 0 & 0 & 1
\end{bmatrix} 
```

## Optical Tracking
Optical tracking consists of two stages where we repeatedly apply the transformations described above

### Define the coordinate systems
We define landmarks on bones that roughly represent their functional axes. Then attach trackers to each bone, which will be captured by the camera.
1. Digitisation of landmarks
2. Definition of bone coordinate system
3. Definition of tracker coordinate system
4. Calculation of bone-to-tracker transform

### Process tracked data
Apply the motion (translations and rotations) captured through the trackers to the landmarks. Describe the motion of the bones relative to each other.
1. Import tracked data
2. Define dynamic tracker positions
3. Calculate dynamic bone positions

The definition of the coordinate systems comes from `create_landmarks()`, and calculation of the transforms is in `bone_to_tracker_transform()`.
## Naming convention
Variable names convey the frame of reference, whether it's a transform or a position vector, the bone, and whether it is a datum or a landmark digitisation.

A variable named `gTt0` is trying to convey the following notation:
```math
_gT_{t_0}
```
which means that we define a **transform** ($T$) in the **global** ($_g$) frame of reference that describes the **tibia**'s ($t$) digitisation ($_0$).
In other words, the matrix that describes transforming from **global** to **tibial** frame of reference. 

Similarly, the variable `tt_r_tc` is:
```math
_{\text{tt}}r_{t_c}
```
which means that in the **Tibial Tracker**'s ($_{\text{tt}}$) frame of reference, we define a position vector ($r$) for the **tibia**'s origin ($_t$).
I.e., the fourth column of the matrix `ttTt`, tibial transform in the tibial tracker's frame of reference.

## Visualising a frame of reference change
`tt_T_t = gT_tt\gTt`, the change of frame of reference can be visualised like so:
```math
\begin{align*}
_{tt}T_{t} &= {\left[ _gT_{tt} \right]}^{-1} \cdot \left[_gT_{t}\right] \\
    & = \left[ _{tt}T_{g} \right] \cdot \left[_gT_{t}\right] \\
    & = \left[ _{tt}T_{\cancel{g}} \right] \cdot \left[\cancel{_g}T_{t}\right] \\
    & = \left[ _{tt}T_{t}\right] \\
\end{align*}
```
