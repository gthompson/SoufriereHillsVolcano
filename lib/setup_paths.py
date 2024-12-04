#!/usr/bin/env python

import os

paths = {}
paths['PROJECT_DIR'] = os.path.join(os.getenv('HOME'), 'Developer','MESS2024_Glenn')
#paths['PROJECT_DIR'] = os.path.join(os.path.expanduser('~'), 'work','PROJECTS','MESS2024_Glenn')
paths['CONTINUOUS_DATA_DIR'] = os.path.join(paths['PROJECT_DIR'], 'data', 'continuous')
paths['SDS_DIR'] = os.path.join(paths['CONTINUOUS_DATA_DIR'], 'SDS')
paths['SDS_VEL_DIR'] = os.path.join(paths['CONTINUOUS_DATA_DIR'], 'SDS_VEL')
paths['SDS_DISP_DIR'] = os.path.join(paths['CONTINUOUS_DATA_DIR'], 'SDS_DISP')
paths['SAM_DIR'] = os.path.join(paths['CONTINUOUS_DATA_DIR'], 'SAM')
paths['SAMBINARY_DIR'] = os.path.join(paths['CONTINUOUS_DATA_DIR'], 'SAM', 'binary')
paths['RESPONSE_DIR'] = os.path.join(paths['PROJECT_DIR'], 'data', 'responses')
paths['EVENTS_DIR'] = os.path.join(paths['PROJECT_DIR'], 'data', 'events')
paths['DB_DIR'] = os.path.join(paths['PROJECT_DIR'], 'db')

for key in paths:
    thisdir = paths[key]
    if not os.path.isdir(thisdir):
        print(f"Making directory {thisdir}")
        os.makedirs(thisdir)
    else:
        #print(f"{thisdir} exists")
        pass

if __name__ == "__main__":
    print(paths)
