#!/usr/bin/env python

import makeplot_cru

def main():

      title = ['Surface air temperature']
      vabsname = [['Tair', 'Tair', 't2m', 'tmp']]  
      filename = ['Tair']                          

      for i in range(0, len(title)):  

          makeplot_cru.main(title[i], filename[i], vabsname[i])

#     makeplot.main(a,b,c,d)
#     a,b is any string that used for title and outputname
#     c,d must be consistent with the variable name in RASM or ERA data


if __name__ == "__main__":
        main()
