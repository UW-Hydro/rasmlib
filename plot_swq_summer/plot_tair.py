#!/usr/bin/env python

import makeplot

def main():

#      makeplot.main('Surface air temperature', 'Tair', 'Tair', 't2m')
#      makeplot.main('Skin temperature', 'Radt', 'Radt', 'skt')
      makeplot.main('Snow water equivalent', 'Swq', 'Swq', 'sd')
#     makeplot.main(a,b,c,d)
#     a,b is any string that used for title and outputname
#     c,d must be consistent with the variable name in RASM or ERA data


if __name__ == "__main__":
        main()
