import os
import re
from typing import Union, Tuple, Optional
import numpy as np
import pandas as pd 

def parse_imrotbl(value):
    # Parse the Imec Readout Table (imRo)
    entries = re.findall(r'\((.*?)\)', value)
    
    # Parse header
    probe_type, num_channels = map(int, entries[0].split(','))

    channels = []
    for entry in entries[1:]:
        channel_data = tuple(map(int, entry.split()))
        
        if probe_type in [0, 1020, 1030, 1100, 1120, 1121, 1122, 1123, 1200, 1300]:  # NP 1.0-like
            if len(channel_data) == 6:
                channels.append({
                    'channel': channel_data[0],  # Channel ID
                    'bank': channel_data[1],     # Bank number of the connected electrode
                    'refid': channel_data[2],    # Reference ID index (0=ext, 1=tip, [2..4]=on-shnk-ref)
                    'apgain': channel_data[3],   # AP band gain
                    'lfgain': channel_data[4],   # LF band gain
                    'apfilt': channel_data[5]    # AP hipass filter applied (1=ON)
                })
                # Note: On-shank ref electrodes are {192,576,960}
        elif probe_type in [21, 2003, 2004]:  # NP 2.0, single multiplexed shank
            if len(channel_data) == 4:
                channels.append({
                    'channel': channel_data[0],    # Channel ID
                    'bank_mask': channel_data[1],  # Bank mask (logical OR of {1=bnk-0, 2=bnk-1, 4=bnk-2, 8=bnk-3})
                    'refid': channel_data[2],      # Reference ID index
                    'electrode': channel_data[3],  # Electrode ID (range [0,1279])
                    'apgain': 80
                })
                # Note for Type-21: Reference ID values are {0=ext, 1=tip, [2..5]=on-shnk-ref}
                # On-shank ref electrodes are {127,507,887,1251}
                # Note for Type-2003,2004: Reference ID values are {0=ext, 1=gnd, 2=tip}
                # On-shank reference electrodes are removed from commercial 2B probes
        elif probe_type in [24, 2013, 2014]:  # NP 2.0, 4-shank
            if len(channel_data) == 5:
                channels.append({
                    'channel': channel_data[0],   # Channel ID
                    'shank': channel_data[1],     # Shank ID (with tips pointing down, shank-0 is left-most)
                    'bank': channel_data[2],      # Bank ID
                    'refid': channel_data[3],     # Reference ID index
                    'electrode': channel_data[4], # Electrode ID (range [0,1279] on each shank)
                    'apgain': 80
                })
            # Note for Type-24: Reference ID values are {0=ext, [1..4]=tip[0..3], [5..8]=on-shnk-0, [9..12]=on-shnk-1, [13..16]=on-shnk-2, [17..20]=on-shnk-3}
            # On-shank ref electrodes of any shank are {127,511,895,1279}
            # Note for Type-2013,2014: Reference ID values are {0=ext, 1=gnd, [2..5]=tip[0..3]}
            # On-shank reference electrodes are removed from commercial 2B probes

    return {'probe_type': probe_type, 'num_channels': num_channels, 'channels': channels}


def parse_snschanmap(value):
    # Parse header
    header_match = re.match(r'\((\d+),(\d+),(\d+),(\d+),(\d+)\)', value)
    imec_header_match = re.match(r'\((\d+),(\d+),(\d+)\)', value)
    
    if header_match:
        header = {
            'MN_channels': int(header_match.group(1)),
            'MA_channels': int(header_match.group(2)),
            'mux_channels': int(header_match.group(3)),
            'XA_channels': int(header_match.group(4)),
            'XD_words': int(header_match.group(5))
        }
    elif imec_header_match:
        header = {
            'AP_channels': int(imec_header_match.group(1)),
            'LF_channels': int(imec_header_match.group(2)),
            'SY_channels': int(imec_header_match.group(3))
        }
    else:
        raise ValueError("Invalid header format in snsChanMap")
    
    # Parse channel map
    channel_map = []
    channel_pattern = r'\(([A-Za-z]+\d+);(\d+):(\d+)\)'
    for match in re.finditer(channel_pattern, value):
        channel_map.append({
            'name': match.group(1),
            'channel': int(match.group(2)),
            'order': int(match.group(3))
        })
    
    return {
        'header': header,
        'channel_map': channel_map
    }


def parse_snsgeommap(value):
    """
    Parse the snsGeomMap value for imec probes.
    
    The GeomMap describes how electrodes are arranged on the probe.
    It consists of a header and electrode entries.
    
    Header format: (part_number,shank_count,shank_spacing,per_shank_width)
    Electrode entry format: (s:x:z:u) where:
        s: zero-based shank number (left-most when tips point down)
        x: x-coordinate (um) of electrode center
        z: z-coordinate (um) of electrode center
        u: 0/1 flag indicating if the electrode is "used"
    
    Note: (X,Z) coordinates are relative to each shank's own origin.
    X-origin is the left edge of the shank, Z-origin is the center of the bottom-most electrode row.
    """
    pattern = r'\(([^)]+)\)'
    matches = re.findall(pattern, value)
    
    if not matches:
        raise ValueError("Invalid snsGeomMap format")
    
    # Parse header
    header = matches[0].split(',')
    
    # Parse electrode entries
    geom_data = []
    for entry in matches[1:]:
        s, x, z, u = map(int, entry.split(':'))
        geom_data.append({
            'shank': s,
            'x': x,
            'z': z,
            'used': bool(u)
        })
    
    return {
        'header': {
            'part_number': header[0],
            'shank_count': int(header[1]),
            'shank_spacing': int(header[2]),
            'per_shank_width': int(header[3])
        },
        'electrodes': geom_data
    }


def get_gain(meta):
    """
    Get the gain of the recording.

    Parameters:
    -----------
    meta : dict
        Metadata dictionary containing recording information.

    Returns:
    --------
    numpy.ndarray
        Array of gain values for each channel.

    Raises:
    -------
    ValueError
        If the recording type is not recognized or if required metadata is missing.
    """
    if meta['typeThis'] == 'imec':
        if 'imroTbl' not in meta or 'channels' not in meta['imroTbl']:
            raise ValueError("Missing 'imroTbl' or 'channels' in metadata for imec recording")
        
        if meta['fileName'].endswith('.ap.bin'):
            return np.array([channel['apgain'] for channel in meta['imroTbl']['channels']])
        elif meta['fileName'].endswith('.lf.bin'):
            return np.array([channel['lfgain'] for channel in meta['imroTbl']['channels']])
        else:
            raise ValueError(f"Unrecognized file type for imec recording: {meta['filename']}")
    
    elif meta['typeThis'] == 'nidq':
        if 'snsMnMaXaDw' not in meta or 'niMNGain' not in meta or 'niMAGain' not in meta:
            raise ValueError("Missing required metadata for nidq recording")
        
        n_mn = meta['snsMnMaXaDw']['MN']
        n_ma = meta['snsMnMaXaDw']['MA']
        n_xa = meta['snsMnMaXaDw']['XA']
        n_dw = meta['snsMnMaXaDw']['DW']
        
        gains = np.ones(n_mn + n_ma + n_xa + n_dw)
        gains[:n_mn] = meta['niMNGain']
        gains[n_mn:n_mn + n_ma] = meta['niMAGain']
        return gains
    
    elif meta['typeThis'] == 'obx':
        if 'snsXaDwSy' not in meta:
            raise ValueError("Missing 'snsXaDwSy' in metadata for obx recording")
        n_xa = meta['snsXaDwSy']['XA']
        n_dw = meta['snsXaDwSy']['DW']
        n_sy = meta['snsXaDwSy']['SY']
        gains = np.ones(n_xa + n_dw + n_sy)
        return gains
    else:
        raise ValueError(f"Unrecognized recording type: {meta['typeThis']}")


def get_uV_per_bit(meta):
    AiRangeMax = meta.get('imAiRangeMax') or meta.get('niAiRangeMax') or meta.get('obAiRangeMax')
    MaxInt = meta.get('imMaxInt') or meta.get('niMaxInt') or meta.get('obMaxInt') or 512
    gains = get_gain(meta)
    return 1000000 * AiRangeMax / MaxInt / gains


def get_channel_idx(meta, analog=True):
    """
    Get the channel index slice based on the recording type and whether analog or digital channels are desired.

    Parameters:
    -----------
    meta : dict
        Metadata dictionary containing recording information.
    analog : bool, optional
        If True, return index for analog channels. If False, return index for digital channels. Default is True.

    Returns:
    --------
    slice
        A slice object representing the channel index range.

    Raises:
    -------
    ValueError
        If the recording type is not recognized.
    """
    if meta['typeThis'] == 'imec':
        n_ap = meta['snsApLfSy']['AP']
        n_lf = meta['snsApLfSy']['LF']
        n_sy = meta['snsApLfSy']['SY']
        if analog:
            channel_idx = slice(0, n_ap + n_lf)
        else:
            channel_idx = slice(n_ap + n_lf, n_ap + n_lf + n_sy)
    elif meta['typeThis'] == 'nidq':
        n_mn = meta['snsMnMaXaDw']['MN']
        n_ma = meta['snsMnMaXaDw']['MA']
        n_xa = meta['snsMnMaXaDw']['XA']
        n_dw = meta['snsMnMaXaDw']['DW']
        if analog:
            channel_idx = slice(0, n_mn + n_ma + n_xa)
        else:
            channel_idx = slice(n_mn + n_ma + n_xa, n_mn + n_ma + n_xa + n_dw)
    elif meta['typeThis'] == 'obx':
        n_xa = meta['snsXaDwSy']['XA']
        n_dw = meta['snsXaDwSy']['DW']
        n_sy = meta['snsXaDwSy']['SY']
        if analog:
            channel_idx = slice(0, n_xa)
        else:
            channel_idx = slice(n_xa, n_xa + n_dw + n_sy)
    else:
        raise ValueError(f"Unrecognized recording type: {meta['typeThis']}")
    
    return channel_idx

def read_meta(filename):
    """
    Read SpikeGLX meta data from a file.

    Parameters
    ----------
    filename : str
        Path to the meta file.

    Returns
    -------
    dict
        Dictionary containing the meta data with format:
        variable_name=int / float / str / list / tuple
    """
    filename_meta = filename.replace('.bin', '.meta')
    if not os.path.exists(filename_meta):
        return None
    meta = {}
    with open(filename_meta, 'r') as f:
        for line in f:
            key, value = line.strip().split('=', 1)
            if key.startswith('~'):
                key = key[1:]  # Remove leading '~'
            
            # Parse complex structures
            if key == 'imroTbl':
                meta[key] = parse_imrotbl(value)
            elif key == 'snsChanMap':
                meta[key] = parse_snschanmap(value)
            elif key == 'snsGeomMap':
                meta[key] = parse_snsgeommap(value)
            elif key in ['acqMnMaXaDw', 'snsMnMaXaDw']:
                values = list(map(int, value.split(',')))
                meta[key] = {
                    'MN': values[0],  # Multiplexed Neural
                    'MA': values[1],  # Multiplexed Auxiliary
                    'XA': values[2],  # Auxiliary Analog
                    'DW': values[3]   # Digital Word
                }
            elif key == 'snsApLfSy':
                values = list(map(int, value.split(',')))
                meta[key] = {
                    'AP': values[0],
                    'LF': values[1],
                    'SY': values[2]
                }
            elif key in ['acqXaDwSy', 'snsXaDwSy']:
                values = list(map(int, value.split(',')))
                meta[key] = {
                    'XA': values[0],
                    'DW': values[1],
                    'SY': values[2]
                }
            else:
                # Try to convert to int or float
                try:
                    meta[key] = int(value)
                except ValueError:
                    try:
                        meta[key] = float(value)
                    except ValueError:
                        # If not int or float, keep as string
                        meta[key] = value
    meta['metaFileName'] = filename_meta
    return meta

def get_channel_idx(meta, analog=True):
    """
    Get the channel index slice based on the recording type and whether analog or digital channels are desired.

    Parameters:
    -----------
    meta : dict
        Metadata dictionary containing recording information.
    analog : bool, optional
        If True, return index for analog channels. If False, return index for digital channels. Default is True.

    Returns:
    --------
    slice
        A slice object representing the channel index range.

    Raises:
    -------
    ValueError
        If the recording type is not recognized.
    """
    if meta['typeThis'] == 'imec':
        n_ap = meta['snsApLfSy']['AP']
        n_lf = meta['snsApLfSy']['LF']
        n_sy = meta['snsApLfSy']['SY']
        if analog:
            channel_idx = slice(0, n_ap + n_lf)
        else:
            channel_idx = slice(n_ap + n_lf, n_ap + n_lf + n_sy)
    elif meta['typeThis'] == 'nidq':
        n_mn = meta['snsMnMaXaDw']['MN']
        n_ma = meta['snsMnMaXaDw']['MA']
        n_xa = meta['snsMnMaXaDw']['XA']
        n_dw = meta['snsMnMaXaDw']['DW']
        if analog:
            channel_idx = slice(0, n_mn + n_ma + n_xa)
        else:
            channel_idx = slice(n_mn + n_ma + n_xa, n_mn + n_ma + n_xa + n_dw)
    elif meta['typeThis'] == 'obx':
        n_xa = meta['snsXaDwSy']['XA']
        n_dw = meta['snsXaDwSy']['DW']
        n_sy = meta['snsXaDwSy']['SY']
        if analog:
            channel_idx = slice(0, n_xa)
        else:
            channel_idx = slice(n_xa, n_xa + n_dw + n_sy)
    else:
        raise ValueError(f"Unrecognized recording type: {meta['typeThis']}")
    
    return channel_idx

def get_probe(meta: dict) -> dict:
    """
    Create a dictionary to track probe information for Kilosort4.

    Parameters
    ----------
    meta : dict
        Dictionary containing metadata information.

    Returns
    -------
    dict
        Dictionary with the following keys, all corresponding to NumPy ndarrays:
        'chanMap': the channel indices that are included in the data.
        'xc':      the x-coordinates (in micrometers) of the probe contact centers.
        'yc':      the y-coordinates (in micrometers) of the probe contact centers.
        'kcoords': shank or channel group of each contact (not used yet, set all to 0).
        'n_chan':  the number of channels.
    """
    geom_data = meta['snsGeomMap']
    electrodes = pd.DataFrame(geom_data['electrodes'])
    shank_spacing = geom_data['header']['shank_spacing']
    
    x = electrodes['x'].values.astype(np.float32) + electrodes['shank'].values.astype(np.float32) * shank_spacing
    y = electrodes['z'].values.astype(np.float32)
    connected = electrodes['used'].values.astype(np.bool_)
    
    channel_map = np.array([i['channel'] for i in meta['snsChanMap']['channel_map']], dtype=np.int32)
    
    channel_idx = get_channel_idx(meta)

    return {
        'chanMap': channel_map[channel_idx] + 1,
        'xc': x,
        'yc': y,
        'kcoords': electrodes['shank'].values.astype(np.int32),
        'connected': connected,}