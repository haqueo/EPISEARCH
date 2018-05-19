def find_indices(N,m):
    # first find the semi triangle
    # number in first row is N-2
    # ALL INDICES ARE ARRAY STYLE, 0,1,..,M-1
    # N is how many actual SNPs the search is going through, so it goes (0,1,2),...,(0,1,N-1)
    # m is how many triples each thread will be searching.
    
    m_orig = m
    triangle_cumsum = 0

    current_count = 0
    found = False
    
    for i in range(N-2,0,-1): # i = 4,3,2,1
        
        old_triangle_cumsum = triangle_cumsum
        triangle_cumsum = triangle_cumsum + i*(i+1)/2
        triangle_index = (N-2)-i
        triangle_rows = i
        if (triangle_cumsum >= m):
            
            # now find exact row of this triangle
            # how many times do i need to keep subtracting rows until i've gone too far
            current_count = old_triangle_cumsum
            for j in reversed(range(triangle_rows)): # 3,2,1,0
                    current_count = current_count + (j+1)
                    if(current_count >=m ):
                        #okay. it is the j+1th row
                        # now find how far along this row it is
                        for k in range(j+1):#
                            if (k +1 +  current_count -(j+1)== m):
                                # we're done. Onto the next one.
                                # it's the (triangle_index)'th triangle,
                                # j'th row from the bottom
                                # and the k'th one
                                # so (triangle_index,N-2-j,N-2-j + 1 + k)
                                
                                

                                m = m + m_orig
                                
        
    return 0


if __name__ == "__main__":
    N = 5000
    m = 200000
    

    find_indices(N,m)
