import simtk.openmm as mm
import os

def make_path(filename):
    try:
        path = os.path.split(filename)[0]
        os.makedirs(path)
    except OSError:
        pass

def get_platform():
    try:
        platform = mm.Platform.getPlatformByName("CUDA")
        key = "CUDA_VISIBLE_DEVICES"
        if key in os.environ:
            device_index = os.environ[key]
        else:
            device_index = 0
        platform.setPropertyDefaultValue("CudaDeviceIndex", "%d" % device_index)
    except:
        platform = mm.Platform.getPlatformByName("CPU")

    return platform
