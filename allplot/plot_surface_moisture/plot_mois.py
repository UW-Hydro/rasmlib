#!/usr/bin/env python

import makeplot_mois

def main():

      a = ['Runoff','Precipitation','Evapotranspiration']
      c = ['Runoff','Precipitation','Evap']
      d = ['ro','tp','e']
      b = c

      for i in range(0,3):

          makeplot_mois.main(a[i], b[i], c[i], d[i])

#     makeplot.main(a,b,c,d)
#     a,b is any string that used for title and outputname
#     c,d must be consistent with the variable name in RASM or ERA data


if __name__ == "__main__":
        main()
