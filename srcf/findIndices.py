def find_indices(N,m):
    # first find the semi triangle
    # number in first row is N-2
    # ALL INDICES ARE ARRAY STYLE, 0,1,..,M-1
    # N is how many actual SNPs the search is going through, so it goes (0,1,2),...,(0,1,N-1)
    # m is how many triples each thread will be searching.
    
    m_orig = m
    triangle_cumsum = 0
    current_count = 0
    
    file = open('indexes.txt','w')
    
    newIndex = (0,1,2)
    
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
                                oldIndex = newIndex
                                newIndex = (triangle_index, N-2-j,N-2-j+1+k)
                                
                                file.write(str(oldIndex[0]) +  ',' + str(oldIndex[1]) + ',' + str(oldIndex[2]) +
                                           ',' + str(newIndex[0]) + ',' + str(newIndex[1]) + ',' + str(newIndex[2]) + '\n')
                                m = m + m_orig
    if newIndex != (N-3,N-2,N-1):
        file.write(str(newIndex[0]) +  ',' + str(newIndex[1]) + ',' + str(newIndex[2]) +
                   ',' + str(N-3) + ',' + str(N-2) + ',' + str(N-1) + '\n')
    
                                        
    file.close()
    return 0

if __name__ == "__main__":
    N = 100
    m = 30
    
    find_indices(N,m)
