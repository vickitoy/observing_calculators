#!/usr/bin/env python

import cgi
#import cgitb
#cgitb.enable()

import sys
sys.path.insert(0, '/pathtofiles')
import limitmag

form = cgi.FieldStorage()

if form.getvalue('mode'):
   mode = form.getvalue('mode')
else:
   mode = "Not set"

if mode == 'L':
    v1 = form.getvalue('L_exptime', 0)
    v2 = form.getvalue('L_sn', 0)
    v3 = form.getvalue('L_filt')
    v4 = form.getvalue('L_mode')
    
    if v3 == 'YJ': v3 = 'Y'
    if v3 == 'HK': v3 = 'H'
    
    labels = ['Exposure time (s)', 'Signal-to-noise ratio', 'Filter', 'Observing mode']
    outtext = 'Calculated AB limiting magnitude:'
    
    result,eta = limitmag.limit('h2rg', v3, int(v4), float(v1), float(v2))
    #try:
    #    # limit(detector, band, R, exptime, sn)
    #    result,eta = limitmag.limit('h2rg', v3, int(v4), float(v1), float(v2))
    #except:
    #    print "Content-Type: text/html"
    #    print "<html><body>Error</body></html>
    #    return
    
elif mode == 'E':
    v1 = '0'
    v2 = '0'

    v1 = form.getvalue('E_mag', 0)
    v2 = form.getvalue('E_sn', 0)
    v3 = form.getvalue('E_filt')
    v4 = form.getvalue('E_mode')

    if v3 == 'YJ': v3 = 'Y'
    if v3 == 'HK': v3 = 'H'

    labels = ['AB magnitude:', 'Signal-to-noise ratio', 'Filter', 'Observing mode']
    outtext = 'Calculated exposure time (sec):'


    result,eta = limitmag.exptimecalc('h2rg', v3, int(v4), float(v1),float(v2))
    
    #try:
    #    # exptimecalc(detector, band, R, objmag, signos)
    #    result,eta = limitmag.exptimecalc('h2rg', v3, int(v4), float(v1),float(v2))
    #except:
    #    except:
    #    print "Content-Type: text/html"
    #    print "<html><body>Error</body></html>
    #    return
    
else:
    v1 = 'N/A'
    v2 = 'N/A'
    v3 = 'N/A'
    v4 = 'N/A'  

omodedict = {'0': 'Photometry', '30': 'R = 30', '4000': 'R = 4000'}
modedict = {'L': 'Limiting Magnitude Calculator', 'E': 'Exposure Time Calculator', 
    'Not Set': 'Error'}


print "Content-Type: text/html"     # HTML is following
print                               # blank line, end of headers
print '<html>'
print '<head>'
print '<title>RIMAS %s</title>' % (modedict[mode])
print """
</head>
<body>
"""
print '<h2> %s %.1f</h2>' % (outtext, result)

print '<h3>%s Inputs:</h3>' % (modedict[mode])
print """
<div>
 <table style="margin-left:20px;" width="300px" border="1" align="left" cellpadding="5">
 """
print '<tr align="left">'
print ' <th>%s:</th>' % (labels[0])
print ' <td>%s</td>' % v1
print '</tr>'
print '<tr align="left">'
print ' <th>%s:</th>' % (labels[1])
print ' <td>%s</td>' % v2
print '</tr>'
print '<tr align="left">'
print ' <th>%s:</th>' % (labels[2])
print ' <td>%s</td>' % v3
print '</tr>'
print '<tr align="left">'
print ' <th>%s:</th>' % (labels[3])
print ' <td>%s</td>' % omodedict[v4]
print '</tr>'
print """
</table>
</div>
""" 


print """
<div style = "clear:both;"></div>
<div>
<h3>Assumptions:</h3>
 <table style="margin-left:20px;" width="300px" border="1" align="left" cellpadding="5">
  <tr align="left">
    <th>Telescope diameter:</th>
    <td>428 cm</td>
  </tr>
  <tr align="left">
    <th>F/#:</th>
    <td>2.47</td>
  </tr>
  <tr align="left">
    <th>Pixel size:</th>
    <td>18 um</td>
  </tr>
  <tr align="left">
    <th>Plate scale:</th>
    <td>0.35 "/pix</td>
  </tr>
  <tr align="left">
    <th>Secondary vignetting:</th>
    <td>0.93</td>
  </tr>
  <tr align="left">
    <th>Binning:</th>
    <td>1x1</td>
  </tr>
  <tr align="left">
    <th>Dark current:</th>
    <td>0.005 e-/s</td>
  </tr>
  <tr align="left">
    <th>Read noise:</th>
    <td>6 e-</td>
  </tr>
  <tr align="left">
    <th>Object spread over N pixels:</th>
    <td>4</td>
  </tr>
  <tr align="left">
    <th>Observing mode efficiency:</th>
"""
print  '<td> %s </td>' % (str(eta))
print """  
  </tr>  
</table>
</div>
""" 

print '</body>'
print '</html>'







