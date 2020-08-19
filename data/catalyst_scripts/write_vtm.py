from paraview.simple import *

from paraview import coprocessing

# This sets the name of the output files. The timestep will automatically be
# appended to the name. For example, if output_name was "data", then resultant
# filenames will be data_1.vtm, data_2.vtm, etc. Also, directories will be
# created with names data_1/, data_2/, etc.
output_name = 'simOutagewww'

# The number of padded zeroes to use for output filenames
padding = 4

# What type of compression to use. If it is set to '', then no compression is
# used. The two available compression schemes are zlib and lz4. The zlib scheme
# will compress better, but will take more time. The lz4 scheme will not
# compress as well, but will be faster to compress. In general, the zlib option
# will not add that much more time, and is recommended.
#compression = ''      # no compression
compression = 'ZLib'   # zlib compression
#compression = 'LZ4'   # lz4 compression

# the frequency to output everything
outputfrequency = 1

#--------------------------------------------------------------
# Code generated from cpstate.py to create the CoProcessor.


# ----------------------- CoProcessor definition -----------------------

def CreateCoProcessor():
  def _CreatePipeline(coprocessor, datadescription):
    class Pipeline:
      adaptorinput = coprocessor.CreateProducer( datadescription, "input" )
      grid = adaptorinput.GetClientSideObject().GetOutputDataObject(0)
      filename = None
      if  grid.IsA('vtkImageData') or grid.IsA('vtkUniformGrid'):
        writer_class = servermanager.writers.XMLPImageDataWriter
        filename = output_name + '_%t.pvti'
      elif  grid.IsA('vtkRectilinearGrid'):
        writer_class = servermanager.writers.XMLPRectilinearGridWriter
        filename = output_name + '_%t.pvtr'
      elif  grid.IsA('vtkStructuredGrid'):
        writer_class = servermanager.writers.XMLPStructuredGridWriter
        filename = output_name + '_%t.pvts'
      elif  grid.IsA('vtkPolyData'):
        writer_class = servermanager.writers.XMLPPolyDataWriter
        filename = output_name + '_%t.pvtp'
      elif  grid.IsA('vtkUnstructuredGrid'):
        #aggregateFilterOut = AggregateDataset(Input = adaptorinput)
        #aggregateFilterOut.NumberOfTargetProcesses = 2
        writer_class = servermanager.writers.XMLPUnstructuredGridWriter
        filename = output_name + '_%t.pvtu'
      elif  grid.IsA('vtkUniformGridAMR'):
        writer_class = servermanager.writers.XMLHierarchicalBoxDataWriter
        filename = output_name + '_%t.vthb'
      elif  grid.IsA('vtkMultiBlockDataSet'):
        writer_class = servermanager.writers.XMLMultiBlockDataWriter
        filename = output_name + '_%t.vtm'
      else:
        print("COPROCESSOR: Don't know how to create a writer for a ", grid.GetClassName())


      

      if compression == '':
        #if  grid.IsA('vtkUnstructuredGrid'):
        #writer = writer_class(Input=aggregateFilterOut)
        #else:
        writer = writer_class(Input=adaptorinput)
      else:
        writer = writer_class(Input=adaptorinput, CompressorType=compression)

      if filename:
        #coprocessor.RegisterWriter(writer, filename, freq=outputfrequency)          # For PV 5.2
        coprocessor.RegisterWriter(writer, filename, freq=outputfrequency, paddingamount=padding)  # For PV 5.3 +

    return Pipeline()

  class CoProcessor(coprocessing.CoProcessor):
    def CreatePipeline(self, datadescription):
      self.Pipeline = _CreatePipeline(self, datadescription)

  #print 'COPROCESSOR: BEFORE BUILD OF COPROCESSOR'
  coprocessor = CoProcessor()
  freqs = {'input': [outputfrequency]}
  coprocessor.SetUpdateFrequencies(freqs)
  return coprocessor

#--------------------------------------------------------------
# Global variables that will hold the pipeline for each timestep
# Creating the CoProcessor object, doesn't actually create the ParaView pipeline.
# It will be automatically setup when coprocessor.UpdateProducers() is called the
# first time.
coprocessor = CreateCoProcessor()

#--------------------------------------------------------------
# Enable Live-Visualizaton with ParaView
coprocessor.EnableLiveVisualization(True)


# ---------------------- Data Selection method ----------------------

def RequestDataDescription(datadescription):
    "Callback to populate the request for current timestep"
    global coprocessor
    if datadescription.GetForceOutput() == True:
        # We are just going to request all fields and meshes from the simulation
        # code/adaptor.
        for i in range(datadescription.GetNumberOfInputDescriptions()):
            datadescription.GetInputDescription(i).AllFieldsOn()
            datadescription.GetInputDescription(i).GenerateMeshOn()
        return

    # setup requests for all inputs based on the requirements of the
    # pipeline.
    coprocessor.LoadRequestedData(datadescription)

# ------------------------ Processing method ------------------------

def DoCoProcessing(datadescription):
    "Callback to do co-processing for current timestep"
    global coprocessor

    # Update the coprocessor by providing it the newly generated simulation data
    # If the pipeline hasn't been setup yet, this will setup the pipeline.
    #print 'COPROCESSOR: UPDATE PRODUCERS'
    coprocessor.UpdateProducers(datadescription)

    # Write output data, if appropriate.
    #print 'COPROCESSOR: WRITE VTM FILES'
    coprocessor.WriteData(datadescription);

    # Write image capture (Last arg: rescale lookup table), if appropriate.
    coprocessor.WriteImages(datadescription, rescale_lookuptable=False)

    # Live Visualization, if enabled.
    coprocessor.DoLiveVisualization(datadescription, "localhost", 22222)
