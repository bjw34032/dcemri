mniLR <- readNIfTI(system.file("nifti/avg152T1_LR_nifti.nii.gz", package="dcemriS4"))
par(bg="black")
image(mniLR)
orthographic(mniLR)
