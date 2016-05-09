'''
Created on 09 May 2016

@author: maxhartig

This is a rudimentary Stokes Flow Solver on a collocated grid using the fipy package
'''

from fipy import *
#from fipy.meshes.grid2D import Grid2D

# declaring some parameters

LX=2.0
LY=1.0
NX=100
NY=50
dX= LX/NX
dY= LY/NY
visc = 1
U = 1
pressureRelax = 0.8
velRelax = 0.5
if __name__ == '__main__':
    sweeps = 300
else:
    sweeps = 5
    
# Building of the mesh

mesh = Grid2D(nx=NX, ny=NY, dx=dX, dy=dY)

#declaring the variables 
pressure= CellVariable(mesh=mesh, name='pressure')
pressureCorrection = CellVariable(mesh=mesh)
xVelocity = CellVariable(mesh=mesh, name='X velocity')
yVelocity = CellVariable(mesh=mesh, name='Y velocity')

# assign velocity as a face variable for the Rie-Chow pressure correction:

velocity = FaceVariable(mesh=mesh, rank=1)

# Stokes Equations in the cell centers:

xVelocityEq = DiffusionTerm(coeff=visc) - pressure.grad.dot([1.,0.])
yVelocityEq = DiffusionTerm(coeff=visc) - pressure.grad.dot([0.,1.])


# going for the pressure correction:

ap= CellVariable(mesh=mesh, value=1.)
coeff= 1./ ap.arithmeticFaceValue*mesh._faceAreas * mesh._cellDistances
pressureCorrectionEq = DiffusionTerm(coeff=coeff) - velocity.divergence

# this type of pressure correction would cause oscillations on a collocated grid. So we use Rie-Chow 

from fipy.variables.faceGradVariable import _FaceGradVariable
volume = CellVariable(mesh=mesh, value=mesh.cellVolumes, name='Volume')
contrvolume= volume.arithmeticFaceValue

xVelocity.constrain(0., mesh.facesRight | mesh.facesLeft | mesh.facesBottom)
xVelocity.constrain(U,mesh.facesTop)
yVelocity.constrain(0.,mesh.exteriorFaces)
X,Y = mesh.faceCenters
pressureCorrection.constrain(0., mesh.facesLeft & (Y < dY))

#the viewers:

if __name__=='__main__':
    viewer = Viewer(vars=(pressure,xVelocity,yVelocity,velocity),xmin=0., xmax=2., ymin=0., ymax=1., colorbar=True)
    
for sweep in range(sweeps):
    # solving the Stokes equations for the starred values:
    xVelocityEq.cacheMatrix()
    xres = xVelocityEq.sweep(var=xVelocity, underRelaxation=velRelax)
    xmat = xVelocityEq.matrix
    
    yres = yVelocityEq.sweep(var=yVelocity, underRelaxation=velRelax)
    
    ## update the ap coefficient from the matrix diagonal
    ap[:] = -xmat.takeDiagonal()
    
    # updating the face velocities with the Rhie-Chow correction
    presgrad = pressure.grad # cell pressure gradient
    fpresgrad = _FaceGradVariable(pressure)
    
    velocity[0] = xVelocity.arithmeticFaceValue \
        + contrvolume / ap.arithmeticFaceValue * \
        (presgrad[0].arithmeticFaceValue-fpresgrad[0])
    velocity[1] = yVelocity.arithmeticFaceValue \
        + contrvolume / ap.arithmeticFaceValue * \
        (presgrad[1].arithmeticFaceValue - fpresgrad[1])
    velocity[..., mesh.exteriorFaces.value] = 0.
    velocity[0, mesh.facesTop.value] = U 
    
    # solving of the pressure correction equation:
    pressureCorrectionEq.cacheRHSvector()
    # no correction at left bottom point because must remain at pressure 0
    pres = pressureCorrectionEq.sweep(var=pressureCorrection)
    rhs = pressureCorrectionEq.RHSvector
    
    # updating the pressure with the corrected values:
    pressure.setValue(pressure+pressureRelax*pressureCorrection)
    #updating velocity with the corrected pressure:
    xVelocity.setValue(xVelocity - pressureCorrection.grad[0] / \
                                            ap * mesh.cellVolumes)
    yVelocity.setValue(yVelocity - pressureCorrection.grad[1] / \
                                            ap * mesh.cellVolumes)
    
    if __name__=='__main__':
        if sweep%10 == 0:
            print 'sweep:', sweep, ', x residual:', xres, \
                                    ', y residual:', yres, \
                                    ', p residual:', pres, \
                                    ', continuity: ', max(abs(rhs)) 
            viewer.plot()
            
print numerix.allclose(pressure.globalValue[...,-1], 162.790867927)