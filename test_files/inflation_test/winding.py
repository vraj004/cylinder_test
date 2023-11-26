## TOOLBOX TO CONVERT WINDING SCHEMES

def toVTK(n, elemClass):
    """
        elemClass = 0 for Hexahedral
        elemClass = 1 for Simplex
    """
    elem = ""
    elem_meshio = ""
    if elemClass == 0:
        if n == 8:
            # Iron Numbering 
            #   *  z = 0       x  z = 1         
            #   *  2------3    x  6------7
            #   *  |      |    x  |      |
            #   *  |      |    x  |      |
            #   *  |      |    x  |      |
            #   *  0------1    x  4------5

            # Xi Numbering 
            #   *  z = 0           z = 1    
            #   *  0------4    *  1------5 
            #   *  |      |    *  |      |    
            #   *  |      |    *  |      |  
            #   *  |      |    *  |      |    
            #   *  2------6    *  3------7   

            # VTK Numbering
            #   *  z = 0       x  z = 1         
            #   *  3------2    x  7------6
            #   *  |      |    x  |      |
            #   *  |      |    x  |      |
            #   *  |      |    x  |      |
            #   *  0------1    x  4------5
            elem = "vtkQuad"
            elem_meshio = "hexahedron"
            conversion = [
                2, 6, 4, 0,
                3, 7, 1, 5
            ]
            return elem, elem_meshio, conversion
        # if n == 20:
        #     elem = "vtkQuadraticHexahedron"
        #     elem_meshio = "hexahedron20"
        if n == 27:
            # Iron Numbering 
            #   *  z = 0           z = 0.5         z = 1          
            #   *  6--7 --8     * 15--16--17    x 24--25--26
            #   *  |      |     *  |      |     x  |      |
            #   *  3  4   5     * 12  13  14    x 21  22  23
            #   *  |      |     *  |      |     x  |      |
            #   *  0--1 --2     *  9--10--11    x 18--19--20

            # Xi Numbering 
            #   *  z = 0           z = 0.5         z = 1    
            #   *  0-- 9--18    *  1--10--19    *  2--11--20     
            #   *  |      |     *  |      |     *  |      |       
            #   *  3  12  21    *  4  13  22    *  5  14  23     
            #   *  |      |     *  |      |     *  |      |     
            #   *  6--15--24    *  7--16--25    *  8--17--26      

            # VTK Numbering 
            #   *  z = 0           z = 0.5         z = 1    
            #   *  3--10--2     * 19--23--18    x  7--14--6
            #   *  |      |     *  |      |     x  |      |
            #   * 11  24  9     * 20  26  21    x 15  25  13
            #   *  |      |     *  |      |     x  |      |
            #   *  0-- 8--1     * 16--22--17    x  4--12--5 
            elem = "vtkTriQuadraticHexahedron"
            elem_meshio = "hexahderon27"
            conversion = [
                6, 24, 18, 0, 8, 26, 20, 2, 
                15, 21, 9, 3, 17, 23, 11, 5, 
                7, 25, 19, 1, 4, 22, 16, 10, 
                12, 14, 13
            ]
        return elem, elem_meshio, conversion
    if elemClass == 1:
        return 
        # if n == 4:
        #     elem = "vtkTetra"
        #     elem_meshio = "tetra"
        # if n == 10:
        #     elem = "vtkQuadTetra"
        #     elem_meshio = "tetra10"

