#!/usr/bin/env python

import makeplot_energy

def main():

      a = ['Latent heat','Sensible heat','Net shortwave radiation','Net longwave radiation',\
          'Downward shortwave radiation','downward longwave radiation']
      c = ['Latht','Senht','Swnet','Lwnet','Swin','Lwin']
      d = ['slhf','sshf','ssr','str','ssrd','strd']
      b = c

      for i in range(0,6):  

          makeplot_energy.main(a[i], b[i], c[i], d[i])

#     makeplot.main(a,b,c,d)
#     a,b is any string that used for title and outputname
#     c,d must be consistent with the variable name in RASM or ERA data


if __name__ == "__main__":
        main()
