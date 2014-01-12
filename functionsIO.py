def writeGrid(gridName, init_condition, Nx1, Ny1, Nx2, Ny2, dx1, dy1, dx2, dy2):
  
  Nx = Nx1 + Nx2
  Ny = Ny1 + Ny2
  check_plane = Nx * Ny

  # coordingates of the left bottom corner cell
  startx = - Nx2 * dx2 - Nx1 * dx1 + 0.5 * dx1
  starty = - 0.5 * Ny2 * dy2 - 0.5 * Ny1 * dy1 + 0.5 * dy1

  # compute sawtooth

  # start point of the sawtooth in the cylinder area
  x0 = - Nx2 * dx2 + 0.5 * dx2
  y0 = - 0.5 * Ny2 * dy2 + 0.5 * dy2
  #print x0
  #print y0
  # array to store saw_position
  saw_position = []
  i = 0

  saw_position.append(Nx-1)

  for j in range(1,Ny2-1):
    while (x0 + i * dx2) ** 2 + (y0 + j * dy2) ** 2 > (0.5 * Ny2 * dy2) ** 2:
      i = i + 1
    saw_position.append(i + Nx1)
    #print i + Nx1
  saw_position.append(Nx-1)

  # compute number of vertical interfaces
  v_num = (Nx + 1) * Ny1 
  for i in range(0,Ny2):
    v_num = v_num + saw_position[i] + 1
  #print v_num

  # compute number of horizontal interfaces in the saw area, not include the two Nx
  h_num = 0
  for i in range(0, Ny2 - 1):
    if saw_position[i] > saw_position[i+1]:
      h_num = h_num + saw_position[i]
    else:
      h_num = h_num + saw_position[i+1]
  # print h_num

  # compute number of cells in the saw area
  c_num = 0
  for i in range(0,Ny2):
    c_num = c_num + saw_position[i]

  # compute number of vertical interfaces in the saw area
  vs_num = c_num + Ny2

  # print startx
  # print starty

  f = open("".join([gridName,'.txt']),"w")

  #write base parameters

  # for num in fluid_params:
  #   f.write(''.join([str(num),'\n']))

  # number of cells
  f.write(''.join([str(check_plane), '\n']))

  # write the cells
  f.write('cells \n')
  for j in range(0, Ny):
    for i in range(0, Nx):
      if i < Nx1 and j < Ny1/2:
        # first 4 are interface id: left, right, bottom, up
        left = j * (Nx + 1) + i
        right = left + 1
        bottom = v_num + j * Nx + i
        top = v_num + (j + 1) * Nx + i
        x = startx + i * dx1
        y = starty + j * dy1
        cellid = j * Nx + i
        mylist = [left, right, bottom, top, x, y, dx1, dy1, cellid]
        for a in range(0,9):
            f.write(str(mylist[a]))
            f.write(', ')
        for a in range(0,3):
            f.write(str(init_condition[a]))
            f.write(',')
        f.write(str(init_condition[3]))
        f.write('\n')
      elif i >= Nx1 and j < Ny1/2:
        left = j * (Nx + 1) + i
        right = left + 1
        bottom = v_num + j * Nx + i
        top = v_num + (j + 1) * Nx + i
        x = startx + (Nx1 - 1) * dx1 + (dx1 + dx2)/2 + (i - Nx1) * dx2
        y = starty + j * dy1
        cellid = j * Nx + i
        mylist = [left, right, bottom, top, x, y, dx2, dy1, cellid]
        for a in range(0,9):
            f.write(str(mylist[a]))
            f.write(', ')
        for a in range(0,3):
            f.write(str(init_condition[a]))
            f.write(',')
        f.write(str(init_condition[3]))
        f.write('\n')
      elif i < Nx1 and j >= Ny1/2 and j < Ny2 + Ny1/2: 
        if j == Ny1/2:
          left = Ny1/2 * (Nx + 1) + i
          right = left + 1
          bottom = v_num + Ny1/2 * Nx + i
          top = bottom + Nx
          x = startx + i * dx1
          y = starty + (Ny1/2 - 1) * dy1 + (dy1 + dy2)/2
          cellid = Ny1/2 * Nx + i
          #print cellid
          mylist = [left, right, bottom, top, x, y, dx1, dy2, cellid]
          for a in range(0,9):
              f.write(str(mylist[a]))
              f.write(', ')
          for a in range(0,3):
              f.write(str(init_condition[a]))
              f.write(',')
          f.write(str(init_condition[3]))
          f.write('\n')
        else:
          temp = 0
          for p in range(0,j - Ny1/2):
            temp = temp + saw_position[p]
          h1_num = 0
          for p in range(0,(j - Ny1/2)):
            if saw_position[p] > saw_position[p+1]:
              h1_num = h1_num + saw_position[p]
            else: 
              h1_num = h1_num + saw_position[p+1]
          left = Ny1/2 * (Nx + 1) + temp + (j-Ny1/2) + i
          right = left + 1
          if saw_position[(j - Ny1/2)-1] > saw_position[(j - Ny1/2)]:
            bottom = v_num + Ny1/2 * Nx + Nx + h1_num - saw_position[(j - Ny1/2)-1] + i
          else:
            bottom = v_num + Ny1/2 * Nx + Nx + h1_num - saw_position[(j - Ny1/2)] + i
          top = v_num + Ny1/2 * Nx + Nx + h1_num + i
          x = startx + i * dx1
          y = starty + (Ny1/2 - 1) * dy1 + (dy1 + dy2)/2 + (j - Ny1/2) * dy2
          cellid = Ny1/2 * Nx + temp + i
          #print cellid
          mylist = [left, right, bottom, top, x, y, dx1, dy2, cellid]
          for a in range(0,9):
              f.write(str(mylist[a]))
              f.write(', ')
          for a in range(0,3):
              f.write(str(init_condition[a]))
              f.write(',')
          f.write(str(init_condition[3]))
          f.write('\n')
      # in the cylinder 
      elif i >= Nx1 and j >= Ny1/2 and j < Ny2 + Ny1/2:
        if saw_position[j - Ny1/2] < i + 1:
          pass
        else:
          if saw_position[j - Ny1/2] - Nx1 > 0:
            # the first row
            if j == Ny1/2:
              for k in range(0,saw_position[0] - Nx1):
                left = Ny1/2 * (Nx + 1) + Nx1 + k
                right = left + 1
                bottom = v_num + Ny1/2 * Nx + Nx1 + k
                top = v_num + Ny1/2 * Nx + Nx + Nx1 + k
                x = x0 + k * dx2
                y = y0 + (j - Ny1/2) * dy2
                cellid = Ny1/2 * Nx + Nx1 + k
                #print cellid
                mylist = [left, right, bottom, top, x, y, dx2, dy2, cellid] 
                for a in range(0,9):
                  f.write(str(mylist[a]))
                  f.write(', ')
                for a in range(0,3):
                  f.write(str(init_condition[a]))
                  f.write(',')
                f.write(str(init_condition[3]))
                f.write('\n')
            # the rest rows
            else:
              #compute horizontal interface in the saw area passed
              h1_num = 0
              for p in range(0,(j - Ny1/2)):
                if saw_position[p] > saw_position[p+1]:
                  h1_num = h1_num + saw_position[p]
                else: 
                  h1_num = h1_num + saw_position[p+1]
              for k in range(0,saw_position[(j - Ny1/2)] - Nx1):
                temp = 0
                for p in range(0,j - Ny1/2):
                  temp = temp + saw_position[p]
                left = Ny1/2 * (Nx + 1) + Nx1 + k + temp + (j - Ny1/2)
                right = left + 1
                if saw_position[(j - Ny1/2)-1] > saw_position[(j - Ny1/2)]:
                  bottom = v_num + Ny1/2 * Nx + Nx + h1_num - saw_position[(j - Ny1/2)-1] + Nx1 + k
                else:
                  bottom = v_num + Ny1/2 * Nx + Nx + h1_num - saw_position[(j - Ny1/2)] + Nx1 + k
                top =  v_num + Ny1/2 * Nx + Nx + h1_num + Nx1 + k
                x = x0 + k * dx2
                y = y0 + (j - Ny1/2) * dy2  
                cellid = Ny1/2 * Nx + Nx1 + k + temp
                #print cellid
                mylist = [left, right, bottom, top, x, y, dx2, dy2, cellid]
                for a in range(0,9):
                  f.write(str(mylist[a]))
                  f.write(', ')
                for a in range(0,3):
                  f.write(str(init_condition[a]))
                  f.write(',')
                f.write(str(init_condition[3]))
                f.write('\n') 
          else:
            pass          
      elif i < Nx1 and j >= Ny2 + Ny1/2:
        temp = 0
        for p in range(0,Ny2):
          temp = temp + saw_position[p]
        #print temp
        left = Ny1/2 * (Nx + 1) + temp + Ny2 + (j - Ny1/2 - Ny2) * (Nx + 1) + i
        right = left + 1
        bottom = v_num + (Ny1/2 + 1) * Nx + h_num + (j - Ny1/2 - Ny2) * Nx + i
        top = bottom + Nx
        x = startx + i * dx1
        y = starty + (Ny1/2 - 1) * dy1 + (dy1 + dy2)/2 + (Ny2 - 1) * dy2 + (dy2 + dy1)/2 + (j - Ny1/2 - Ny2) * dy1 
        cellid = Ny1/2 * Nx + temp + (j - Ny1/2 - Ny2) * Nx + i
        #print cellid
        mylist = [left, right, bottom, top, x, y, dx1, dy1, cellid]
        for a in range(0,9):
            f.write(str(mylist[a]))
            f.write(', ')
        for a in range(0,3):
            f.write(str(init_condition[a]))
            f.write(',')
        f.write(str(init_condition[3]))
        f.write('\n')
      else:
        temp = 0
        for p in range(0,Ny2):
          temp = temp + saw_position[p]
        left = Ny1/2 * (Nx + 1) + temp + Ny2 + (j - Ny1/2 - Ny2) * (Nx + 1) + i
        right = left + 1
        bottom = v_num + Ny1/2 * Nx + Nx + h_num + (j - Ny1/2 - Ny2) * Nx + i
        #print j
        #print bottom
        top = bottom + Nx
        x = startx + (Nx1 - 1) * dx1 + 0.5 * (dx1 + dx2) + (i - Nx1) * dx2
        y = starty + (Ny1/2 - 1) * dy1 + (dy1 + dy2)/2 + (Ny2 - 1) * dy2 + (dy2 + dy1)/2 + (j - Ny1/2 - Ny2) * dy1  
        cellid = Ny1/2 * Nx + temp + (j - Ny1/2 - Ny2) * Nx + i
        mylist = [left, right, bottom, top, x, y, dx2, dy1, cellid]
        for a in range(0,9):
            f.write(str(mylist[a]))
            f.write(', ')
        for a in range(0,3):
            f.write(str(init_condition[a]))
            f.write(',')
        f.write(str(init_condition[3]))
        f.write('\n')

  # write the interfaces
  f.write('interfaces \n')

  # vertical interfaces
  for j in range(0, Ny):
    for i in range(0, Nx + 1):
      if j < Ny1/2:
        if i == 0:
          cellid1 = -1
          cellid2 = j * Nx
          interfaceid = j * (Nx + 1) + i
          mylist = [cellid1, cellid2, interfaceid]
          for a in range(0,2):
            f.write(str(mylist[a]))
            f.write(',')
          f.write(str(mylist[2]))
          f.write(str('\n'))
        elif i == Nx:
          cellid1 = j * Nx + Nx - 1
          cellid2 = -1
          interfaceid = j * (Nx + 1) + i
          mylist = [cellid1, cellid2, interfaceid]
          for a in range(0,2):
            f.write(str(mylist[a]))
            f.write(',')
          f.write(str(mylist[2]))
          f.write(str('\n'))
        else:
          cellid1 = j * Nx + i - 1
          cellid2 = j * Nx + i
          interfaceid = j * (Nx + 1) + i
          mylist = [cellid1, cellid2, interfaceid]
          for a in range(0,2):
            f.write(str(mylist[a]))
            f.write(',')
          f.write(str(mylist[2]))
          f.write(str('\n'))
      elif j >= Ny1/2 and j < Ny1/2 + Ny2:
        if i > saw_position[j-Ny1/2]:
          pass
        else:
          temp = 0
          for p in range(0, j-Ny1/2 + 1):
            temp = temp + saw_position[p]
          if i == 0:
            cellid1 = -1
            cellid2 = Ny1/2 * Nx + temp - saw_position[j-Ny1/2]
            interfaceid = Ny1/2 * (Nx + 1) + temp - saw_position[j-Ny1/2] + j - Ny1/2
            mylist = [cellid1, cellid2, interfaceid]
            for a in range(0,2):
              f.write(str(mylist[a]))
              f.write(',')
            f.write(str(mylist[2]))
            f.write(str('\n'))
          elif i == saw_position[j-Ny1/2]:
            cellid1 = Ny1/2 * Nx + temp - 1
            cellid2 = -3
            interfaceid = Ny1/2 * (Nx + 1) + temp - saw_position[j-Ny1/2] + j - Ny1/2 + saw_position[j-Ny1/2]
            mylist = [cellid1, cellid2, interfaceid]
            for a in range(0,2):
              f.write(str(mylist[a]))
              f.write(',')
            f.write(str(mylist[2]))
            f.write(str('\n'))
          else:
            cellid1 = Ny1/2 * Nx + temp - saw_position[j-Ny1/2] + i - 1
            cellid2 = Ny1/2 * Nx + temp - saw_position[j-Ny1/2] + i
            interfaceid = Ny1/2 * (Nx + 1) + temp - saw_position[j-Ny1/2] + j - Ny1/2 + i
            mylist = [cellid1, cellid2, interfaceid]
            for a in range(0,2):
              f.write(str(mylist[a]))
              f.write(',')
            f.write(str(mylist[2]))
            f.write(str('\n'))
      # the region above the saw
      else:
        # number of vertical interfaces in the middle region
        v1_num = 0
        for p in range(0,Ny2):
          v1_num = v1_num + saw_position[p] + 1
        if i == 0:
          cellid1 = -1
          cellid2 = Ny1/2 * Nx + v1_num - Ny2 + (j - Ny1/2 - Ny2) * Nx
          interfaceid = Ny1/2 * (Nx + 1) + v1_num + (j - Ny1/2 - Ny2) * (Nx + 1)
          mylist = [cellid1, cellid2, interfaceid]
          for a in range(0,2):
            f.write(str(mylist[a]))
            f.write(',')
          f.write(str(mylist[2]))
          f.write(str('\n'))
        elif i == Nx:
          cellid1 = Ny1/2 * Nx + v1_num - Ny2 + (j - Ny1/2 - Ny2) * Nx + Nx - 1
          cellid2 = -1
          interfaceid = Ny1/2 * (Nx + 1) + v1_num + (j - Ny1/2 - Ny2) * (Nx + 1) + Nx
          mylist = [cellid1, cellid2, interfaceid]
          for a in range(0,2):
            f.write(str(mylist[a]))
            f.write(',')
          f.write(str(mylist[2]))
          f.write(str('\n'))
        else:
          cellid1 = Ny1/2 * Nx + v1_num - Ny2 + (j - Ny1/2 - Ny2) * Nx + i - 1
          cellid2 = Ny1/2 * Nx + v1_num - Ny2 + (j - Ny1/2 - Ny2) * Nx + i
          interfaceid = Ny1/2 * (Nx + 1) + v1_num + (j - Ny1/2 - Ny2) * (Nx + 1) + i
          mylist = [cellid1, cellid2, interfaceid]
          for a in range(0,2):
            f.write(str(mylist[a]))
            f.write(',')
          f.write(str(mylist[2]))
          f.write(str('\n'))

  # horizontal interfaces
  for j in range(0, Ny + 1):
    for i in range(0, Nx):
      if j < Ny1/2:
        if j == 0:
          cellid1 = -2
          cellid2 = i
          interfaceid = v_num + i
          mylist = [cellid1, cellid2, interfaceid]
          for a in range(0,2):
            f.write(str(mylist[a]))
            f.write(',')
          f.write(str(mylist[2]))
          f.write(str('\n'))
        else:
          cellid1 = (j-1) * Nx + i
          cellid2 = cellid1 + Nx
          interfaceid = v_num + j * Nx + i
          mylist = [cellid1, cellid2, interfaceid]
          for a in range(0,2):
            f.write(str(mylist[a]))
            f.write(',')
          f.write(str(mylist[2]))
          f.write(str('\n'))
      elif j >= Ny1/2 and j <= Ny1/2 + Ny2:
        if j == Ny1/2:
          cellid1 = (j-1) * Nx + i
          if i < saw_position[0]:
            cellid2 = cellid1 + Nx
          else:
            cellid2 = -3
          interfaceid = v_num + Ny1/2 * Nx + i
          mylist = [cellid1, cellid2, interfaceid]
          for a in range(0,2):
            f.write(str(mylist[a]))
            f.write(',')
          f.write(str(mylist[2]))
          f.write(str('\n'))
        elif j == Ny1/2 + Ny2:
          if i < saw_position[Ny2 - 1]:
            cellid1 = Ny1/2 * Nx + c_num - saw_position[Ny2-1] + i
          else:
            cellid1 = -3
          cellid2 = Ny1/2 * Nx + c_num + i
          interfaceid = v_num + (Ny1/2 + 1) * Nx + h_num + i
          #print interfaceid
          mylist = [cellid1, cellid2, interfaceid]
          for a in range(0,2):
            f.write(str(mylist[a]))
            f.write(',')
          f.write(str(mylist[2]))
          f.write(str('\n'))
        else:
          if i > saw_position[j-Ny1/2] - 1 and i > saw_position[j-Ny1/2 - 1] - 1:
            pass
          else:
            # number of cells passed in the saw
            c1_num = 0
            for p in range(0,j-Ny1/2):
              c1_num = c1_num + saw_position[p]
            h1_num = 0
            for p in range(0,j - Ny1/2):
              if saw_position[p] > saw_position[p+1]:
                h1_num = h1_num + saw_position[p]
              else: 
                h1_num = h1_num + saw_position[p+1]
            if saw_position[j-Ny1/2] > saw_position[j-Ny1/2 - 1]:
              if i < saw_position[j-Ny1/2-1]:
                cellid1 = Ny1/2 * Nx + c1_num - saw_position[j-Ny1/2 - 1] + i
              else:
                cellid1 = -3
              cellid2 = Ny1/2 * Nx + c1_num + i
              interfaceid = v_num + (Ny1/2 + 1) * Nx + h1_num - saw_position[j-Ny1/2] + i
              mylist = [cellid1, cellid2, interfaceid]
              for a in range(0,2):
                f.write(str(mylist[a]))
                f.write(',')
              f.write(str(mylist[2]))
              f.write(str('\n'))
            else:
              cellid1 = Ny1/2 * Nx + c1_num - saw_position[j-Ny1/2 - 1] + i
              if i < saw_position[j-Ny1/2]:
                cellid2 = Ny1/2 * Nx + c1_num + i
              else:
                cellid2 = -3
              interfaceid = v_num + (Ny1/2 + 1) * Nx + h1_num - saw_position[j-Ny1/2-1] + i
              mylist = [cellid1, cellid2, interfaceid]
              for a in range(0,2):
                f.write(str(mylist[a]))
                f.write(',')
              f.write(str(mylist[2]))
              f.write(str('\n'))
      else:
        if j < Ny:
          cellid1 = Ny1/2 * Nx + c_num + (j - Ny1/2 - Ny2 - 1) * Nx + i
          cellid2 = cellid1 + Nx
          interfaceid = v_num + (Ny1/2 + 1) * Nx + h_num + (j-Ny1/2-Ny2) * Nx + i
          mylist = [cellid1, cellid2, interfaceid]
          for a in range(0,2):
            f.write(str(mylist[a]))
            f.write(',')
          f.write(str(mylist[2]))
          f.write(str('\n'))
        else:
          cellid1 = Ny1/2 * Nx + c_num + (Ny1/2 - 1) * Nx + i
          cellid2 = -2
          interfaceid = v_num + (Ny1/2 + 1) * Nx + h_num + (j-Ny1/2-Ny2) * Nx + i
          mylist = [cellid1, cellid2, interfaceid]
          for a in range(0,2):
            f.write(str(mylist[a]))
            f.write(',')
          f.write(str(mylist[2]))
          f.write(str('\n'))

  f.close()

  return 0