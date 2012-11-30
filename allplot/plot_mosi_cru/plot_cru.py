#!/usr/bin/env python

import makeplot_cru

def main():

      a = ['Precipitation']
      c = ['Precipitation']
      d = ['tp']
      b = c

      for i in range(0,1):  # now is only for precipitation

          makeplot_cru.main(a[i], b[i], c[i], d[i])

#     makeplot.main(a,b,c,d)
#     a,b is any string that used for title and outputname
#     c,d must be consistent with the variable name in RASM or ERA data


if __name__ == "__main__":
        main()
