############################################################################
#this computes forces on the boundary due to hookes law stiff*displacement
############################################################################


def target_force(b1, b2, bt1, bt2, fb1, fb2, stiff, Q):

    for i in range(Q):
        fb1[i,0] = stiff*(bt1[i,0] - b1[i,0])   #spring constant*displacement between target and boundary points
        fb2[i,0] = stiff*(bt2[i,0] - b2[i,0])

    #print fb1, fb2
    return fb1, fb2 
