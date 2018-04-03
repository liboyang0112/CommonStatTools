import ROOT

def getHistograms(a, outdir=0):
  result = []
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
              h.SetDirectory(outdir)
              result.append(h)
            else:
              print k.GetName(), k.GetClassName() 
              print kk.GetName(), kk.GetClassName() 
              print kkk.GetName(), kkk.GetClassName() 
              raise RuntimeError('Object is not a TH1F or TH1D')
  return result

def run(inputname, outputname, force_overwrite):
  f_in = ROOT.TFile.Open(inputname)
  if force_overwrite:
    f_out = ROOT.TFile.Open(outputname, 'RECREATE')
  else:
    f_out = ROOT.TFile.Open(outputname, 'CREATE')
  
  getHistograms(f_in, f_out)
  f_out.Write()
  f_out.IsA().Destructor(f_out)
  f_in.IsA().Destructor(f_in)

if __name__ == '__main__':
  from argparse import ArgumentParser
  
  parser = ArgumentParser(description='extract histograms from an HistFactory workspace file and put them elsewhere', add_help=True)
  parser.add_argument('-i', '--input', type=str, dest='input', help='path to the workspace file', metavar='input')
  parser.add_argument('-o', '--output', type=str, dest='output', help='path to the output ROOT file', metavar='output')
  parser.add_argument('-f', '--force-overwrite', action='store_true', dest='force_overwrite', help='overwrite the output ROOT file')
  
  options = parser.parse_args()
  
  if (options.input == None):
    parser.error('Must specify old and new files')
  if (options.output == None):
    parser.error('Must specify old and new files')
  
  run(options.input, options.output, options.force_overwrite)

