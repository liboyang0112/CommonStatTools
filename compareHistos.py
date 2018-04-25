import ROOT

def getDirectoryReport(a):
  report = {}
  for k in a.GetListOfKeys():
    if k.GetClassName() == 'TDirectoryFile':
      # found a channel; now loop over samples
      channel = k.ReadObj()
      for kk in channel.GetListOfKeys():
        if kk.GetClassName() == 'TDirectoryFile':
          # found a sample; now loop over histograms
          sample = kk.ReadObj()
          for kkk in sample.GetListOfKeys():
            if kkk.GetClassName() == 'TH1F' or kkk.GetClassName() == 'TH1D':
              h = kkk.ReadObj()
              name = '%s_%s_%s' % (k.GetName(), kk.GetName(), kkk.GetName())
#             name = name.replace('lowMass', 'S2only') # to compare with Xinran
              h.SetDirectory(0)
              report[name] = h
            else:
              print k.GetName(), k.GetClassName() 
              print kk.GetName(), kk.GetClassName() 
              print kkk.GetName(), kkk.GetClassName() 
              raise RuntimeError('Object is not a TH1F or TH1D')
  return report

def compareValues(name, val_a, val_b, tolerance=0):
  if val_b == 0:
    return '%s: %.3f vs %.3f' % (name, val_a, val_b)
  if abs(val_a/val_b-1) > tolerance:
    return '%s: %.3f vs %.3f (%0.2f%%)' % (name, val_a, val_b, (val_b/val_a-1.)*100.)
  else:
    return '%s: identical' % (name)

def compareFiles(a, b):
  rep_a = getDirectoryReport(a)
  rep_b = getDirectoryReport(b)

  only_a = filter(lambda x: x not in rep_b, rep_a.keys())
  only_b = filter(lambda x: x not in rep_a, rep_b.keys())
  common = filter(lambda x: x in rep_a, rep_b.keys())
  
  print '\nElements only in %s:' % a.GetName()
  for el in only_a: print '  %s' % el

  print '\nElements only in %s:' % b.GetName()
  for el in only_b: print '  %s' % el

  print '\nCommon elements:'
  for el in sorted(common):
    printout = [ '  %s' % el ]
    printout.append( compareValues('number of bins', rep_a[el].GetNbinsX(), rep_b[el].GetNbinsX()) )
    printout.append( compareValues('x-axis min', rep_a[el].GetXaxis().GetXmin(), rep_b[el].GetXaxis().GetXmin()) )
    printout.append( compareValues('x-axis max', rep_a[el].GetXaxis().GetXmax(), rep_b[el].GetXaxis().GetXmax()) )
    printout.append( compareValues('integral', rep_a[el].Integral(), rep_b[el].Integral()) )
    printout.append( compareValues('mean', rep_a[el].GetMean(), rep_b[el].GetMean()) )
    printout.append( compareValues('rms', rep_a[el].GetRMS(), rep_b[el].GetRMS()) )

    print '\n    '.join(printout)

      

if __name__ == '__main__':
    from argparse import ArgumentParser

    parser = ArgumentParser(description='compare histograms from two HistFactory workspace files', add_help=True)
    parser.add_argument('-o', '--old', type=str, dest='old', help='path to the old workspace file', metavar='old')
    parser.add_argument('-n', '--new', type=str, dest='new', help='path to the new workspace file', metavar='new')

    options = parser.parse_args()

    if (options.old == None):
      parser.error('Must specify old and new files')
    if (options.new == None):
      parser.error('Must specify old and new files')

    old = ROOT.TFile.Open(options.old)
    new = ROOT.TFile.Open(options.new)

    compareFiles(old, new)
