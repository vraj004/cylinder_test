#A module listing functions to generate simple meshes.

import math
import numpy as np

def annulus(r_inner,r_outer,z_length,numberOfRadialElements,numberOfCircumferentialElements,numberOfZElements,InterpolationType):
    #lists that will be created and returned.
    #set up nodes.
    node_list = []
    node_idx = 0
    node_idx_list = []
    bottom_node_list = []
    top_node_list = []
    yfix_node_list = []
    xfix_node_list = []
    internal_node_list = []
    outer_node_list = []

    #origin
    xorig = 0.0 # x origin
    yorig = 0.0 # y origin
    zorig = 0.0 # z origin
    theta_orig = 2 * math.pi # theta origin

    elem_r_thickness = (r_outer - r_inner)/numberOfRadialElements
    elem_theta_rad = 2.0*math.pi / numberOfCircumferentialElements
    elem_z_height = z_length / numberOfZElements
    if (InterpolationType == 1):
        r_delta = elem_r_thickness
        ridx_end = numberOfRadialElements+1
        theta_delta = elem_theta_rad
        thetaidx_end = numberOfCircumferentialElements
        z_delta = elem_z_height
        zidx_end = numberOfZElements+1
    elif (InterpolationType == 2):
        r_delta = elem_r_thickness/2.0
        ridx_end = numberOfRadialElements+2
        theta_delta = elem_theta_rad / 2.0
        thetaidx_end = numberOfCircumferentialElements*2
        z_delta = elem_z_height/2.0
        zidx_end = numberOfZElements+2

    for ridx in range(0, ridx_end, 1):
        r = r_inner + ridx * r_delta
        print(r)
        for thetaidx in range(0, thetaidx_end, 1):
            theta = theta_orig - thetaidx * theta_delta
            for zidx in range(0, zidx_end, 1):
                z = zorig + zidx * z_delta
                node_list.append(
                    [
                        xorig + r * math.cos(theta),
                        yorig + r * math.sin(theta),
                        zorig + z
                    ]
                )
                node_idx = node_idx + 1 #node numbering starting at 1.
                node_idx_list.append([node_idx])
                print('node_idx=',node_idx)
                print([
                    xorig + r * math.cos(theta),
                    yorig + r * math.sin(theta),
                    zorig + z
                ])
                if (thetaidx == 0 or thetaidx == 4):
                    yfix_node_list.append([node_idx])
                if (zidx == 0):
                    bottom_node_list.append([node_idx])
                if (zidx == zidx_end-1):
                    top_node_list.append([node_idx])
                if (thetaidx == 2 or thetaidx == 6):
                    xfix_node_list.append([node_idx])
                if (ridx == 0):
                    internal_node_list.append([node_idx])
                if (ridx==ridx_end-1):
                    outer_node_list.append([node_idx])
    print('Total nodes in mesh=',len(node_list))
    #set up elements
    elemidx = 0
    elem_start_theta2pi_nodeidx = 0
    elem_start_nodeidx = 0
    for zelemidx in range(0, numberOfZElements):

        for radialelemidx in range(0, numberOfRadialElements):
            for circumelemidx in range(0, numberOfCircumferentialElements):
                elemidx = elemidx + 1  # elem numbering starting at 1.
                if(InterpolationType == 1):
                    print('Linear interpolated element definition not implemented yet')
                if (InterpolationType == 2):
                    print('elem_start_nodeidx=', elem_start_nodeidx)
                    elem_startmid_nodeidx = elem_start_nodeidx + (2 * 3 * numberOfCircumferentialElements)
                    print('elem_startmid_nodeidx=', elem_startmid_nodeidx)
                    elem_startend_nodeidx = elem_startmid_nodeidx + (2 * 3 * numberOfCircumferentialElements)
                    print('elem_startend_nodeidx', elem_startend_nodeidx)
                    if(circumelemidx==numberOfCircumferentialElements-1):
                        elem_startmid_theta2pi_nodeidx = elem_start_theta2pi_nodeidx + (2 * 3 * numberOfCircumferentialElements)
                        elem_startend_theta2pi_nodeidx = elem_startmid_theta2pi_nodeidx + (2 * 3 * numberOfCircumferentialElements)

                        elem_nodesxi3startidx = list(range(elem_start_nodeidx, elem_start_nodeidx + 6, 1))+list(range(elem_start_theta2pi_nodeidx, elem_start_theta2pi_nodeidx + 3, 1))
                        elem_nodesxi3mididx = list(range(elem_startmid_nodeidx, elem_startmid_nodeidx + 6, 1))+list(range(elem_startmid_theta2pi_nodeidx, elem_startmid_theta2pi_nodeidx + 3, 1))
                        elem_nodesxi3endidx = list(range(elem_startend_nodeidx, elem_startend_nodeidx + 6, 1))+list(range(elem_startend_theta2pi_nodeidx, elem_startend_theta2pi_nodeidx + 3, 1))
                    else:
                        elem_nodesxi3startidx = list(range(elem_start_nodeidx, elem_start_nodeidx + 9, 1))
                        elem_nodesxi3mididx = list(range(elem_startmid_nodeidx, elem_startmid_nodeidx + 9, 1))
                        elem_nodesxi3endidx = list(range(elem_startend_nodeidx, elem_startend_nodeidx + 9, 1))

                    elem_node_indices = elem_nodesxi3startidx + elem_nodesxi3mididx + elem_nodesxi3endidx

                    # elem_node_list.append()
                    print('Element number', elemidx)
                    print(elem_nodesxi3startidx)
                    print(elem_nodesxi3mididx)
                    print(elem_nodesxi3endidx)
                    print(elem_node_indices)
                    elem_start_nodeidx = elem_start_nodeidx + 6
    elem_start_theta2pi_nodeidx = elem_start_theta2pi_nodeidx + 6
    return node_list, node_idx_list, top_node_list,bottom_node_list, yfix_node_list, xfix_node_list, internal_node_list, outer_node_list
