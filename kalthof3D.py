import gmsh
import math

def create_v_cut_rectangle(ix, iy, TH, theta,wz,imx,imy,imz, distance_from_corner=0):
    gmsh.initialize()
    gmsh.model.add("rectangle_v_cut_occ")

    # Create the rectangle
    rectangle = gmsh.model.occ.addRectangle(0, 0, 0, ix, iy)

    # Compute the base of the isosceles triangle
    triangle_base = 2 * TH * math.tan(math.radians(theta) / 2)

    # Points and lines for the first triangle (left corner)
    p1 = gmsh.model.occ.addPoint(distance_from_corner, iy, 0)
    p2 = gmsh.model.occ.addPoint(distance_from_corner + triangle_base/2, iy - TH, 0)
    p3 = gmsh.model.occ.addPoint(distance_from_corner + triangle_base, iy, 0)
    l1 = gmsh.model.occ.addLine(p1, p2)
    l2 = gmsh.model.occ.addLine(p2, p3)
    l3 = gmsh.model.occ.addLine(p3, p1)
    triangle1 = gmsh.model.occ.addPlaneSurface([gmsh.model.occ.addWire([l1, l2, l3])])

    # Points and lines for the second triangle (right corner)
    p4 = gmsh.model.occ.addPoint(ix - distance_from_corner, iy, 0)
    p5 = gmsh.model.occ.addPoint(ix - distance_from_corner - triangle_base/2, iy - TH, 0)
    p6 = gmsh.model.occ.addPoint(ix - distance_from_corner - triangle_base, iy, 0)
    l4 = gmsh.model.occ.addLine(p4, p5)
    l5 = gmsh.model.occ.addLine(p5, p6)
    l6 = gmsh.model.occ.addLine(p6, p4)
    triangle2 = gmsh.model.occ.addPlaneSurface([gmsh.model.occ.addWire([l4, l5, l6])])

    # Cut the triangles from the rectangle
    remaining, _ = gmsh.model.occ.cut([(2, rectangle)], [(2, triangle1)])
    gmsh.model.occ.cut(remaining, [(2, triangle2)])
    gmsh.model.occ.extrude([(2,remaining[0][1])],0,0,wz)
    
    # Create impactor on top of the plate
    #impactor = gmsh.model.occ.addBox((ix-imx)/2,(ix-imy)/2,wz,imx,imy,imz)
    dmp = 5
    box_start_x = (ix-imx)/2 
    box_start_y = iy+dmp
    box_start_z = (wz-imz)/2
    impactor = gmsh.model.occ.addBox(box_start_x,box_start_y,box_start_z,imx,imy,imz)

    gmsh.model.occ.synchronize()
    gmsh.model.mesh.generate(3)
    gmsh.write("rectangle_v_cut_occ.brep")
    gmsh.fltk.run()
    gmsh.finalize()
    

# Example usage:
ix, iy = 30, 20
imx,imy,imz = 7,10,7
wz =5
TH =8 
theta =3 
distance_from_corner = 10  # Distance of triangles from the rectangle's corners
create_v_cut_rectangle(ix, iy, TH, theta,wz,imx,imy,imz, distance_from_corner)

