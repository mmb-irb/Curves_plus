      function dotdelta(m)
      include 'curves_data.inc'
      common/bisect/bsx(3)
      common/scr/uint(n3*100,4,3),urot(4,3),usc(6),var(6),theta,dist
      dx=bsx(1)-uint(m,4,1)
      dy=bsx(2)-uint(m,4,2)
      dz=bsx(3)-uint(m,4,3)
      r = sqrt(dx*dx+dy*dy+dz*dz)
      dotdelta=(dx*uint(m,3,1)+dy*uint(m,3,2)+dz*uint(m,3,3))/r
      return
      end
