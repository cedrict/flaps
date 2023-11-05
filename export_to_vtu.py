###############################################################################

import numpy as np
from basis_functions import *
from analytical_solution import *
from gravity_vector import *
from density import *
from viscosity import *

###############################################################################


def export_solution_to_vtu(istep,NV,nel,xV,zV,iconV,u,v,vr,vt,q,vel_unit,rad,theta,nx,nz,sr2,\
                           density_nodal,density_elemental,viscosity_nodal,viscosity_elemental,R1,R2,rho_m,gravity_model,\
                           g0,rhoc,rhoblob,Rblob,zblob,hull,inner_element,outer_element,\
                           innerQ2,outerQ2,bc_fix,exp,e_rr2,e_tt2,e_rt2):

   vtufile=open("solution_"+str(istep)+".vtu","w")
   vtufile.write("<VTKFile type='UnstructuredGrid' version='0.1' byte_order='BigEndian'> \n")
   vtufile.write("<UnstructuredGrid> \n")
   vtufile.write("<Piece NumberOfPoints=' %5d ' NumberOfCells=' %5d '> \n" %(NV,nel))
   #####
   vtufile.write("<Points> \n")
   vtufile.write("<DataArray type='Float32' NumberOfComponents='3' Format='ascii'> \n")
   for i in range(0,NV):
       vtufile.write("%e %e %e \n" %(xV[i],0,zV[i]))
   vtufile.write("</DataArray>\n")
   vtufile.write("</Points> \n")
   #####
   vtufile.write("<PointData Scalars='scalars'>\n")
   #--
   vtufile.write("<DataArray type='Float32' NumberOfComponents='3' Name='velocity (cm/year)' Format='ascii'> \n")
   for i in range(0,NV):
       vtufile.write("%e %e %e \n" %(u[i]/vel_unit,0,v[i]/vel_unit))
   vtufile.write("</DataArray>\n")
   #--
   vtufile.write("<DataArray type='Float32' Name='pressure (q)' Format='ascii'> \n")
   for i in range(0,NV):
       vtufile.write("%e \n" %q[i])
   vtufile.write("</DataArray>\n")
   #--
   vtufile.write("<DataArray type='Float32' Name='r' Format='ascii'> \n")
   for i in range(0,NV):
       vtufile.write("%e \n" %rad[i])
   vtufile.write("</DataArray>\n")
   #--
   vtufile.write("<DataArray type='Float32' Name='theta (co-latitude)' Format='ascii'> \n")
   for i in range(0,NV):
       vtufile.write("%e \n" %theta[i])
   vtufile.write("</DataArray>\n")
   #--
   vtufile.write("<DataArray type='Float32' Name='sr2' Format='ascii'> \n")
   for i in range(0,NV):
       vtufile.write("%e \n" %sr2[i])
   vtufile.write("</DataArray>\n")
   #--
   vtufile.write("<DataArray type='Float32' NumberOfComponents='1' Name='density' Format='ascii'> \n")
   for i in range(0,NV):
       vtufile.write("%e \n" %density_nodal[i])
   vtufile.write("</DataArray>\n")
   #--
   vtufile.write("<DataArray type='Float32' Name='viscosity' Format='ascii'> \n")
   for i in range(0,NV):
       vtufile.write("%e \n" %viscosity_nodal[i])
   vtufile.write("</DataArray>\n")
   #--
   vtufile.write("<DataArray type='Float32' NumberOfComponents='3' Name='gravity' Format='ascii'> \n")
   for i in range(0,NV):
       g_x,g_z=gravity_acceleration(xV[i],zV[i],R1,R2,gravity_model,g0,rho_m,rhoc,rhoblob,Rblob,zblob)
       vtufile.write("%e %e %e \n" %(g_x,0,g_z))
   vtufile.write("</DataArray>\n")
   #--
   vtufile.write("<DataArray type='Float32' NumberOfComponents='3' Name='velocity(r,theta)' Format='ascii'> \n")
   for i in range(0,NV):
       vtufile.write("%e %e %e \n" %(vr[i]/vel_unit,0,vt[i]/vel_unit))
   vtufile.write("</DataArray>\n")
   #--
   vtufile.write("<DataArray type='Float32' NumberOfComponents='3' Name='normal vector' Format='ascii'> \n")
   for i in range(0,NV):       
       vtufile.write("%e %e %e \n" %(nx[i],0,nz[i]))
   vtufile.write("</DataArray>\n")
   #--
   vtufile.write("<DataArray type='Float32' Name='hull' Format='ascii'> \n")
   for i in range(0,NV):
       if hull[i]:
          vtufile.write("%e \n" % 1)
       else:
          vtufile.write("%e \n" % 0)
   vtufile.write("</DataArray>\n")
   #--
   vtufile.write("<DataArray type='Float32' Name='outerQ2' Format='ascii'> \n")
   for i in range(0,NV):
       if outerQ2[i]:
          vtufile.write("%e \n" % 1)
       else:
          vtufile.write("%e \n" % 0)
   vtufile.write("</DataArray>\n")
   #--
   vtufile.write("<DataArray type='Float32' Name='innerQ2' Format='ascii'> \n")
   for i in range(0,NV):
       if innerQ2[i]:
          vtufile.write("%e \n" % 1)
       else:
          vtufile.write("%e \n" % 0)
   vtufile.write("</DataArray>\n")
   #--
   vtufile.write("<DataArray type='Float32' Name='bc_fix(u)' Format='ascii'> \n")
   for i in range(0,NV):
       if bc_fix[2*i]:
          vtufile.write("%e \n" % 1)
       else:
          vtufile.write("%e \n" % 0)
   vtufile.write("</DataArray>\n")
   #--
   vtufile.write("<DataArray type='Float32' Name='bc_fix(v)' Format='ascii'> \n")
   for i in range(0,NV):
       if bc_fix[2*i+1]:
          vtufile.write("%e \n" % 1)
       else:
          vtufile.write("%e \n" % 0)
   vtufile.write("</DataArray>\n")
   #--
   vtufile.write("<DataArray type='Float32' Name='err' Format='ascii'> \n")
   for i in range(0,NV):
       vtufile.write("%e \n" %e_rr2[i])
   vtufile.write("</DataArray>\n")
   #--
   vtufile.write("<DataArray type='Float32' Name='ett' Format='ascii'> \n")
   for i in range(0,NV):
       vtufile.write("%e \n" %e_tt2[i])
   vtufile.write("</DataArray>\n")
   #--
   vtufile.write("<DataArray type='Float32' Name='ert' Format='ascii'> \n")
   for i in range(0,NV):
       vtufile.write("%e \n" %e_rt2[i])
   vtufile.write("</DataArray>\n")
   #--
   if exp==0:
      vtufile.write("<DataArray type='Float32' NumberOfComponents='3' Name='velocity(th)' Format='ascii'> \n")
      for i in range(0,NV):
          vtufile.write("%e %e %e \n" %(velocity_x(xV[i],zV[i],R1,R2,exp),0,velocity_y(xV[i],zV[i],R1,R2,exp)))
      vtufile.write("</DataArray>\n")
      #--
      vtufile.write("<DataArray type='Float32' Name='pressure (th)' Format='ascii'> \n")
      for i in range (0,NV):
          vtufile.write("%e \n" % pressure(xV[i],zV[i],R1,R2,rho_m,g0,exp))
      vtufile.write("</DataArray>\n")
      #--
      vtufile.write("<DataArray type='Float32' Name='exx (th)' Format='ascii'> \n")
      for i in range(0,NV):
          vtufile.write("%e \n" %(sr_xx(xV[i],zV[i],R1,R2,exp)))
      vtufile.write("</DataArray>\n")
      #--
      vtufile.write("<DataArray type='Float32' Name='eyy (th)' Format='ascii'> \n")
      for i in range(0,NV):
          vtufile.write("%e \n" %(sr_yy(xV[i],zV[i],R1,R2,exp)))
      vtufile.write("</DataArray>\n")
      #--
      vtufile.write("<DataArray type='Float32' Name='exy (th)' Format='ascii'> \n")
      for i in range(0,NV):
          vtufile.write("%e \n" %(sr_xy(xV[i],zV[i],R1,R2,exp)))
      vtufile.write("</DataArray>\n")
   #--
   vtufile.write("</PointData>\n")
   #####
   vtufile.write("<CellData Scalars='scalars'>\n")
   #--
   vtufile.write("<DataArray type='Float32' Name='viscosity' Format='ascii'> \n")
   for iel in range(0,nel):
       vtufile.write("%e \n" %viscosity_elemental[iel])
   vtufile.write("</DataArray>\n")
   #--
   vtufile.write("<DataArray type='Float32' Name='density' Format='ascii'> \n")
   for iel in range(0,nel):
       vtufile.write("%e \n" %density_elemental[iel])
   vtufile.write("</DataArray>\n")
   #--
   vtufile.write("<DataArray type='Float32' Name='outer_element' Format='ascii'> \n")
   for iel in range(0,nel):
       if outer_element[iel]:
          vtufile.write("%e \n" % 1)
       else:
          vtufile.write("%e \n" % 0)
   vtufile.write("</DataArray>\n")
   #
   vtufile.write("<DataArray type='Float32' Name='inner_element' Format='ascii'> \n")
   for iel in range(0,nel):
       if inner_element[iel]:
          vtufile.write("%e \n" % 1)
       else:
          vtufile.write("%e \n" % 0)
   vtufile.write("</DataArray>\n")
   #--
   vtufile.write("</CellData>\n")
   #####
   vtufile.write("<Cells>\n")
   #--
   vtufile.write("<DataArray type='Int32' Name='connectivity' Format='ascii'> \n")
   for iel in range (0,nel):
       vtufile.write("%d %d %d %d %d %d %d %d %d\n" %(iconV[0,iel],iconV[1,iel],iconV[2,iel],iconV[3,iel],\
                                                      iconV[4,iel],iconV[5,iel],iconV[6,iel],iconV[7,iel],iconV[8,iel]))
   vtufile.write("</DataArray>\n")
   #--
   vtufile.write("<DataArray type='Int32' Name='offsets' Format='ascii'> \n")
   for iel in range (0,nel):
       vtufile.write("%d \n" %((iel+1)*9))
   vtufile.write("</DataArray>\n")
   #--
   vtufile.write("<DataArray type='Int32' Name='types' Format='ascii'>\n")
   for iel in range (0,nel):
       vtufile.write("%d \n" %28)
   vtufile.write("</DataArray>\n")
   #--
   vtufile.write("</Cells>\n")
   #####
   vtufile.write("</Piece>\n")
   vtufile.write("</UnstructuredGrid>\n")
   vtufile.write("</VTKFile>\n")
   vtufile.close()

   #--
   #if compute_sr1:
   #   vtufile.write("<DataArray type='Float32' Name='sr1' Format='ascii'> \n")
   #   for i in range(0,NV):
   #       vtufile.write("%e \n" %sr1[i])
   #   vtufile.write("</DataArray>\n")
   #if compute_sr3:
   #   vtufile.write("<DataArray type='Float32' Name='sr3' Format='ascii'> \n")
   #   for i in range(0,NV):
   #       vtufile.write("%e \n" %sr3[i])
   #   vtufile.write("</DataArray>\n")

###############################################################################

def export_quadrature_points_to_vtu(nqperdim,nqel,qcoords_r,qcoords_s,mapping,xmapping,ymapping):

   vtufile=open("quadrature_points_"+str(nqperdim)+".vtu","w")
   vtufile.write("<VTKFile type='UnstructuredGrid' version='0.1' byte_order='BigEndian'> \n")
   vtufile.write("<UnstructuredGrid> \n")
   vtufile.write("<Piece NumberOfPoints='%5d' NumberOfCells='%5d'> \n" %(nqel,nqel))
   #####
   vtufile.write("<Points> \n")
   vtufile.write("<DataArray type='Float32' NumberOfComponents='3' Format='ascii'> \n")
   for k in range(0,nqel):
       rq=qcoords_r[k]
       sq=qcoords_s[k]
       NNNV=NNN(rq,sq,mapping)
       xq=np.dot(NNNV[:],xmapping[:,0])
       yq=np.dot(NNNV[:],ymapping[:,0])
       vtufile.write("%e %e %e \n" %(xq,yq,0.))
   vtufile.write("</DataArray>\n")
   vtufile.write("</Points> \n")
   #####
   vtufile.write("<Cells>\n")
   #--
   vtufile.write("<DataArray type='Int32' Name='connectivity' Format='ascii'> \n")
   for iel in range (0,nqel):
       vtufile.write("%d \n" %(iel))
   vtufile.write("</DataArray>\n")
   #--
   vtufile.write("<DataArray type='Int32' Name='offsets' Format='ascii'> \n")
   for iel in range (0,nqel):
       vtufile.write("%d \n" %((iel+1)*1))
   vtufile.write("</DataArray>\n")
   #--
   vtufile.write("<DataArray type='Int32' Name='types' Format='ascii'>\n")
   for iel in range (0,nqel):
       vtufile.write("%d \n" %1)
   vtufile.write("</DataArray>\n")
   #--
   vtufile.write("</Cells>\n")
   #####
   vtufile.write("</Piece>\n")
   vtufile.write("</UnstructuredGrid>\n")
   vtufile.write("</VTKFile>\n")
   vtufile.close()

###############################################################################

def export_mapping_points_to_vtu(mapping,mmapping,xmapping,ymapping):

   vtufile=open("mapping_points_"+mapping+".vtu","w")
   vtufile.write("<VTKFile type='UnstructuredGrid' version='0.1' byte_order='BigEndian'> \n")
   vtufile.write("<UnstructuredGrid> \n")
   vtufile.write("<Piece NumberOfPoints='%5d' NumberOfCells='%5d'> \n" %(mmapping,mmapping))
   #####
   vtufile.write("<Points> \n")
   vtufile.write("<DataArray type='Float32' NumberOfComponents='3' Format='ascii'> \n")
   for i in range(0,mmapping):
       vtufile.write("%e %e %e \n" %(xmapping[i,0],ymapping[i,0],0.))
   vtufile.write("</DataArray>\n")
   vtufile.write("</Points> \n")
   #####
   vtufile.write("<Cells>\n")
   #--
   vtufile.write("<DataArray type='Int32' Name='connectivity' Format='ascii'> \n")
   for iel in range (0,mmapping):
       vtufile.write("%d \n" %(iel))
   vtufile.write("</DataArray>\n")
   #--
   vtufile.write("<DataArray type='Int32' Name='offsets' Format='ascii'> \n")
   for iel in range (0,mmapping):
       vtufile.write("%d \n" %((iel+1)*1))
   vtufile.write("</DataArray>\n")
   #--
   vtufile.write("<DataArray type='Int32' Name='types' Format='ascii'>\n")
   for iel in range (0,mmapping):
       vtufile.write("%d \n" %1)
   vtufile.write("</DataArray>\n")
   #--
   vtufile.write("</Cells>\n")
   #####
   vtufile.write("</Piece>\n")
   vtufile.write("</UnstructuredGrid>\n")
   vtufile.write("</VTKFile>\n")
   vtufile.close()

###############################################################################

def export_Q1_mesh_to_vtu(NV,nel,xV,yV,iconQ1):

   vtufile=open("mesh_Q1.vtu","w")
   vtufile.write("<VTKFile type='UnstructuredGrid' version='0.1' byte_order='BigEndian'> \n")
   vtufile.write("<UnstructuredGrid> \n")
   vtufile.write("<Piece NumberOfPoints=' %5d ' NumberOfCells=' %5d '> \n" %(NV,nel))
   #####
   vtufile.write("<Points> \n")
   vtufile.write("<DataArray type='Float32' NumberOfComponents='3' Format='ascii'> \n")
   for i in range(0,NV):
       vtufile.write("%e %e %e \n" %(xV[i],yV[i],0.))
   vtufile.write("</DataArray>\n")
   vtufile.write("</Points> \n")
   #####
   vtufile.write("<Cells>\n")
   #--
   vtufile.write("<DataArray type='Int32' Name='connectivity' Format='ascii'> \n")
   for iel in range (0,nel):
       vtufile.write("%d %d %d %d\n" %(iconQ1[0,iel],iconQ1[1,iel],iconQ1[2,iel],iconQ1[3,iel]))
   vtufile.write("</DataArray>\n")
   #--
   vtufile.write("<DataArray type='Int32' Name='offsets' Format='ascii'> \n")
   for iel in range (0,nel):
       vtufile.write("%d \n" %((iel+1)*4))
   vtufile.write("</DataArray>\n")
   #--
   vtufile.write("<DataArray type='Int32' Name='types' Format='ascii'>\n")
   for iel in range (0,nel):
       vtufile.write("%d \n" %9)
   vtufile.write("</DataArray>\n")
   #--
   vtufile.write("</Cells>\n")
   #####
   vtufile.write("</Piece>\n")
   vtufile.write("</UnstructuredGrid>\n")
   vtufile.write("</VTKFile>\n")
   vtufile.close()

###############################################################################

def export_gravity_to_vtu(istep,np_grav,xM,zM,gvect_x,gvect_z):

   filename = 'gravity_{:04d}.vtu'.format(istep)
   vtufile=open(filename,"w")
   vtufile.write("<VTKFile type='UnstructuredGrid' version='0.1' byte_order='BigEndian'> \n")
   vtufile.write("<UnstructuredGrid> \n")
   vtufile.write("<Piece NumberOfPoints='%5d' NumberOfCells='%5d'> \n" %(np_grav,np_grav))
   vtufile.write("<Points> \n")
   vtufile.write("<DataArray type='Float32' NumberOfComponents='3' Format='ascii'> \n")
   for i in range(0,np_grav):
       vtufile.write("%e %e %e \n" %(xM[i],zM[i],0.))
   vtufile.write("</DataArray>\n")
   vtufile.write("</Points> \n")
   vtufile.write("<PointData Scalars='scalars'>\n")
   vtufile.write("<DataArray type='Float32' NumberOfComponents='3' Name='gravity' Format='ascii'> \n")
   for i in range(0,np_grav):
       vtufile.write("%e %e %e \n" %(gvect_x[i],gvect_z[i],0))
   vtufile.write("</DataArray>\n")
   vtufile.write("</PointData>\n")
   vtufile.write("<Cells>\n")
   vtufile.write("<DataArray type='Int32' Name='connectivity' Format='ascii'> \n")
   for i in range (0,np_grav):
       vtufile.write("%d \n" %(i))
   vtufile.write("</DataArray>\n")
   vtufile.write("<DataArray type='Int32' Name='offsets' Format='ascii'> \n")
   for i in range (0,np_grav):
       vtufile.write("%d \n" %(i+1))
   vtufile.write("</DataArray>\n")
   vtufile.write("<DataArray type='Int32' Name='types' Format='ascii'>\n")
   for iel in range (0,np_grav):
       vtufile.write("%d \n" %1)
   vtufile.write("</DataArray>\n")
   vtufile.write("</Cells>\n")
   vtufile.write("</Piece>\n")
   vtufile.write("</UnstructuredGrid>\n")
   vtufile.write("</VTKFile>\n")
   vtufile.close()

###############################################################################

def export_slices(nel_phi,NV,nel,r,theta,iconV,rho):
    dphi=2*np.pi/nel_phi
    
    phi=np.zeros(NV,dtype=np.float64) 
    for jel in range(0,nel_phi):

        phi[:]=jel*dphi

        xxx=r[:]*np.sin(theta[:])*np.cos(phi[:])
        yyy=r[:]*np.sin(theta[:])*np.sin(phi[:])
        zzz=r[:]*np.cos(theta[:])

        vtufile=open("slice_"+str(jel)+".vtu","w")
        vtufile.write("<VTKFile type='UnstructuredGrid' version='0.1' byte_order='BigEndian'> \n")
        vtufile.write("<UnstructuredGrid> \n")
        vtufile.write("<Piece NumberOfPoints=' %5d ' NumberOfCells=' %5d '> \n" %(NV,nel))
        #####
        vtufile.write("<Points> \n")
        vtufile.write("<DataArray type='Float32' NumberOfComponents='3' Format='ascii'> \n")
        for i in range(0,NV):
            vtufile.write("%e %e %e \n" %(xxx[i],yyy[i],zzz[i]))
        vtufile.write("</DataArray>\n")
        vtufile.write("</Points> \n")
        #####
        vtufile.write("<CellData Scalars='scalars'>\n")
        vtufile.write("<DataArray type='Float32' NumberOfComponents='1' Name='density' Format='ascii'> \n")
        for iel in range(0,nel):
            vtufile.write("%e \n" % rho[iel])
        vtufile.write("</DataArray>\n")
        vtufile.write("</CellData>\n")
        #####
        vtufile.write("<Cells>\n")
        #--
        vtufile.write("<DataArray type='Int32' Name='connectivity' Format='ascii'> \n")
        for iel in range (0,nel):
            vtufile.write("%d %d %d %d %d %d %d %d\n" %(iconV[0,iel],iconV[1,iel],iconV[2,iel],iconV[3,iel],\
                                                        iconV[4,iel],iconV[5,iel],iconV[6,iel],iconV[7,iel]))
        vtufile.write("</DataArray>\n")
        #--
        vtufile.write("<DataArray type='Int32' Name='offsets' Format='ascii'> \n")
        for iel in range (0,nel):
            vtufile.write("%d \n" %((iel+1)*8))
        vtufile.write("</DataArray>\n")
        #--
        vtufile.write("<DataArray type='Int32' Name='types' Format='ascii'>\n")
        for iel in range (0,nel):
            vtufile.write("%d \n" %23)
        vtufile.write("</DataArray>\n")
        #--
        vtufile.write("</Cells>\n")
        #####
        vtufile.write("</Piece>\n")
        vtufile.write("</UnstructuredGrid>\n")
        vtufile.write("</VTKFile>\n")
        vtufile.close()

###############################################################################
