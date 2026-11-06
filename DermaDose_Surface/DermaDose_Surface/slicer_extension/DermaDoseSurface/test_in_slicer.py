"""Integration test: execute with Slicer --no-main-window --testing --python-script."""
import sys
import time
import traceback
from pathlib import Path
import tempfile
import numpy as np
import qt
import slicer

root = Path(__file__).resolve().parent
sys.path.insert(0, str(root))
import DermaDoseSurface as module


def run():
    test = module.DermaDoseSurfaceTest()
    test.runTest()
    data = np.full((48,64,64), -1000, dtype=np.int16)
    data[3:45,4:60,5:59] = 50
    volume = slicer.util.addVolumeFromArray(data, name='DermaDose async test')
    parent = slicer.qMRMLWidget()
    parent.setLayout(qt.QVBoxLayout())
    widget = module.DermaDoseSurfaceWidget(parent)
    widget.setup()
    widget.volume.setCurrentNode(volume)
    ticks = [0]
    timer = qt.QTimer()
    timer.setInterval(25)
    timer.connect('timeout()', lambda: ticks.__setitem__(0, ticks[0] + 1))
    timer.start()
    widget.onExtract()
    assert widget._worker is not None
    deadline = time.perf_counter() + 120
    while widget._worker is not None and time.perf_counter() < deadline:
        slicer.app.processEvents()
        time.sleep(.01)
    assert widget._worker is None, 'Worker timed out'
    assert widget.logic.segmentationNode is not None, widget.report.plainText
    assert ticks[0] > 2, 'The Qt event loop did not respond'
    with tempfile.TemporaryDirectory() as folder:
        meta = widget.logic.export(Path(folder)/'surface.stl')
        assert meta['stl_quality']['boundary_edges'] == 0
        assert meta['stl_quality']['components'] == 1
    previous = widget.logic.segmentationNode
    widget.onExtract()
    widget.onCancelExtraction()
    deadline = time.perf_counter() + 15
    while widget._worker is not None and time.perf_counter() < deadline:
        slicer.app.processEvents()
        time.sleep(.01)
    assert widget._worker is None
    assert widget.logic.segmentationNode is previous
    timer.stop()
    widget.cleanup()
    print('DermaDose Surface integration PASS; Qt ticks:', ticks[0], flush=True)


if __name__ == '__main__':
    try:
        run()
        slicer.app.exit(0)
    except Exception:
        traceback.print_exc()
        slicer.app.exit(1)
