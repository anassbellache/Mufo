# MUFO

Minimal MATLAB helpers for tomography preprocessing and reconstruction.

## Contents

- `mufo/simple_preprocess.m` – flat-field correction, sinogram generation and
  center-of-rotation detection.
- `mufo/simple_reconstruct.m` – basic filtered backprojection using `iradon`.
- `mufo/tasks/CenterOfRotation.m` – helper class used for center detection.
- `mufo/tests/` – unit tests.

## Quick start

```matlab
proj   = rand(64,64,180);   % raw projections
flat   = ones(64,64);       % flat-field reference
dark   = zeros(64,64);      % dark-field reference
angles = linspace(0, pi, 180);

[corrected, sino, center] = mufo.simple_preprocess(proj, flat, dark, angles);
volume = mufo.simple_reconstruct(sino, angles, center);
```

Run the tests to verify the toolbox:

```matlab
run('run_tests.m')
```

These functions are implemented purely in MATLAB and do not require the
complete UFO installation.

