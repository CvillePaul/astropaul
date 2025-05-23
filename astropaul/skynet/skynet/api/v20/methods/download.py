"""
SkyNet API procedural Python interface: download methods
"""

from __future__ import absolute_import, division, print_function

from ..client import api_call


__all__ = ['jpg', 'png', 'fits', 'header', 'movie']


def jpg(api_key=None, server=None, **kwargs):
    """
    Download image(s) in JPEG format

    :param str api_key: optional API access key
    :param str server: optional API server (for testing)
    :param kwargs::
        image | images: requested image ID(s); multiple IDs are provided as a
            list or a single string, with individual IDs separated by commas.
            An image ID consists of a single-character image type identifier,
            followed by an image ID number. Type Identifiers:
                "r" => raw exposure,
                "m" => master calibration image,
                "w" => Afterglow workspace image,
                "t" => Afterglow temp image,
                "s" => Afterglow sample image
        obs | observations: requested optical observation ID(s); multiple IDs
            are provided as a list or a single string, with individual IDs
            separated by commas
        layer: for multi-HDU observations (e.g. polarimetry or spectral), return
            specific layer (0-based); if omitted, return the primary image
            (layer 0)
        total_parts: split the ZIP archive returned into multiple parts and
            return a single part per request; see also `part`. Ignored for a
            single-image request.
        part: if `total_parts` is set and > 1, this is the number of part to
            return, starting from 1. Ignored for a single-image request.
        wcs = 1: only images with world coordinate system (WCS) in their headers
            will be returned. Ignored for a single-image request.
        filter: when present, only images of the specified filter will be
            returned. Ignored for a single-image request.
        explen: when present, only images of the specified exposure length
            (seconds) will be returned. Ignored for a single-image request.
        reduce: bias, dark, and flat correct the requested image(s)
        reducequiet: same as "reduce" except that in the event of a failure
            during reduction the unreduced image is returned.
        delta: maximum separation (in days) allowed between master calibration
            images and the source image during reduction
        mbias: override the auto-selection of master bias to be used in
            calibration by specifying its ID. Presence of this variable will
            automatically add the option "reduce". If set to 0, this type of
            calibration will not be applied.
        mdark: override the auto-selection of master dark to be used in
            calibration by specifying its ID. Presence of this variable will
            automatically add the option "reduce". If set to 0, this type of
            calibration will not be applied.
        mflat: override the auto-selection of master flat to be used in
            calibration by specifying its ID.  Presence of this variable will
            automatically add the option "reduce". If set to 0, this type of
            calibration will not be applied.
        remove_cosmics = 1: remove cosmic rays from the images; defaults to 1 if
            reduce=1 and to 0 otherwise
        find_cosmics = 1: return images containing only cosmic rays present in
            the original images.
        bin: uniform software binning factor.
        hbin: horizontal software binning factor.
        vbin: vertical software binning factor.
        scale: scale images by the given factor
        width: scale images to the given width
        height: scale images to the given width
        min: specifies the lower percentile of the histogram used to set the
            black point of the image
        max: specifies the upper percentile of the histogram used to set the
            white point of the image
        cmap: produce a false-color image; should be one of the colormap names
            supported by matplotlib; see
            http://matplotlib.org/users/colormaps.html;
            default: "gray" (produce a grayscale image)
        quality: quality of JPEG compression

    :return: a pair (filename, data); for a single-image request, `data`
        contains the JPEG image; for a multiple-image request, `data` is the
        ZIP archive contents
    """
    return api_call(
        'download/jpg', 'get', kwargs, server=server, api_key=api_key)


def png(api_key=None, server=None, **kwargs):
    """
    Download image(s) in PNG format

    :param str api_key: optional API access key
    :param str server: optional API server (for testing)
    :param kwargs::
        image | images: requested image ID(s); multiple IDs are provided as a
            list or a single string, with individual IDs separated by commas.
            An image ID consists of a single-character image type identifier,
            followed by an image ID number. Type Identifiers:
                "r" => raw exposure,
                "m" => master calibration image,
                "w" => Afterglow workspace image,
                "t" => Afterglow temp image,
                "s" => Afterglow sample image
        obs | observations: requested optical observation ID(s); multiple IDs
            are provided as a list or a single string, with individual IDs
            separated by commas
        layer: for multi-HDU observations (e.g. polarimetry or spectral), return
            specific layer (0-based); if omitted, return the primary image
            (layer 0)
        total_parts: split the ZIP archive returned into multiple parts and
            return a single part per request; see also `part`. Ignored for a
            single-image request.
        part: if `total_parts` is set and > 1, this is the number of part to
            return, starting from 1. Ignored for a single-image request.
        wcs = 1: only images with world coordinate system (WCS) in their headers
            will be returned. Ignored for a single-image request.
        filter: when present, only images of the specified filter will be
            returned. Ignored for a single-image request.
        explen: when present, only images of the specified exposure length
            (seconds) will be returned. Ignored for a single-image request.
        reduce: bias, dark, and flat correct the requested image(s)
        reducequiet: same as "reduce" except that in the event of a failure
            during reduction the unreduced image is returned.
        delta: maximum separation (in days) allowed between master calibration
            images and the source image during reduction
        mbias: override the auto-selection of master bias to be used in
            calibration by specifying its ID. Presence of this variable will
            automatically add the option "reduce". If set to 0, this type of
            calibration will not be applied.
        mdark: override the auto-selection of master dark to be used in
            calibration by specifying its ID. Presence of this variable will
            automatically add the option "reduce". If set to 0, this type of
            calibration will not be applied.
        mflat: override the auto-selection of master flat to be used in
            calibration by specifying its ID.  Presence of this variable will
            automatically add the option "reduce". If set to 0, this type of
            calibration will not be applied.
        remove_cosmics = 1: remove cosmic rays from the images; defaults to 1 if
            reduce=1 and to 0 otherwise
        find_cosmics = 1: return images containing only cosmic rays present in
            the original images.
        scale: scale images by the given factor
        width: scale images to the given width
        height: scale images to the given width
        min: specifies the lower percentile of the histogram used to set the
            black point of the image
        max: specifies the upper percentile of the histogram used to set the
            white point of the image
        cmap: produce a false-color image; should be one of the colormap names
            supported by matplotlib; see
            http://matplotlib.org/users/colormaps.html;
            default: "gray" (produce a grayscale image)

    :return: a pair (filename, data); for a single-image request, `data`
        contains the PNG image; for a multiple-image request, `data` is the
        ZIP archive contents
    """
    return api_call(
        'download/png', 'get', kwargs, server=server, api_key=api_key)


def fits(api_key=None, server=None, **kwargs):
    """
    Download image(s) in FITS format

    :param str api_key: optional API access key
    :param str server: optional API server (for testing)
    :param kwargs::
        image | images: requested image ID(s); multiple IDs are provided as a
            list or a single string, with individual IDs separated by commas.
            An image ID consists of a single-character image type identifier,
            followed by an image ID number. Type Identifiers:
                "r" => raw exposure,
                "m" => master calibration image,
                "w" => Afterglow workspace image,
                "t" => Afterglow temp image,
                "s" => Afterglow sample image
        obs | observations: requested optical observation ID(s); multiple IDs
            are provided as a list or a single string, with individual IDs
            separated by commas
        radio_obs: radio observation ID(s); multiple IDs are provided as a list
            or a single string, with individual IDs separated by commas
        layer: for multi-HDU observations (e.g. polarimetry, spectral, or
            radio), return specific layer (0-based); if omitted, return
            the whole FITS file
        total_parts: split the ZIP archive returned into multiple parts and
            return a single part per request; see also `part`. Ignored for a
            single-image request.
        part: if `total_parts` is set and > 1, this is the number of part to
            return, starting from 1. Ignored for a single-image request.
        wcs = 1: only images with world coordinate system (WCS) in their headers
            will be returned. Ignored for a single-image request.
        filter: when present, only images of the specified filter will be
            returned. Ignored for a single-image request.
        explen: when present, only images of the specified exposure length
            (seconds) will be returned. Ignored for a single-image request.
        reduce: bias, dark, and flat correct the requested image(s)
        reducequiet: same as "reduce" except that in the event of a failure
            during reduction the unreduced image is returned.
        delta: maximum separation (in days) allowed between master calibration
            images and the source image during reduction
        mbias: override the auto-selection of master bias to be used in
            calibration by specifying its ID. Presence of this variable will
            automatically add the option "reduce". If set to 0, this type of
            calibration will not be applied.
        mdark: override the auto-selection of master dark to be used in
            calibration by specifying its ID. Presence of this variable will
            automatically add the option "reduce". If set to 0, this type of
            calibration will not be applied.
        mflat: override the auto-selection of master flat to be used in
            calibration by specifying its ID.  Presence of this variable will
            automatically add the option "reduce". If set to 0, this type of
            calibration will not be applied.
        remove_cosmics = 1: remove cosmic rays from the images; defaults to 1 if
            reduce=1 and to 0 otherwise
        find_cosmics = 1: return images containing only cosmic rays present in
            the original images.
        scale: scale images by the given factor
        width: scale images to the given width
        height: scale images to the given width
        force_int: convert FITS images returned to unsigned 16 bit
        processed: return specific radio cartographer channel(s)
        channel: if `processed` is set, return the given radio cartographer
            channel ("left", "right", or "composite"); default: return all three
            channels

    :return: a pair (filename, data); for a single-image request, `data`
        contains the FITS image; for a multiple-image request, `data` is the
        ZIP archive contents
    """
    return api_call(
        'download/fits', 'get', kwargs, server=server, api_key=api_key)


def header(api_key=None, server=None, **kwargs):
    """
    Download image header(s)

    :param str api_key: optional API access key
    :param str server: optional API server (for testing)
    :param kwargs::
        image | images: requested image ID(s); multiple IDs are provided as a
            list or a single string, with individual IDs separated by commas.
            An image ID consists of a single-character image type identifier,
            followed by an image ID number. Type Identifiers:
                "r" => raw exposure,
                "m" => master calibration image,
                "w" => Afterglow workspace image,
                "t" => Afterglow temp image,
                "s" => Afterglow sample image
        obs | observations: requested optical observation ID(s); multiple IDs
            are provided as a list or a single string, with individual IDs
            separated by commas
        radio_obs: radio observation ID(s); multiple IDs are provided as a list
            or a single string, with individual IDs separated by commas
        layer: for multi-HDU observations (e.g. polarimetry or spectral), return
            specific layer (0-based); if omitted, return the primary HDU header
        total_parts: split the ZIP archive returned into multiple parts and
            return a single part per request; see also `part`. Ignored for a
            single-image request.
        part: if `total_parts` is set and > 1, this is the number of part to
            return, starting from 1. Ignored for a single-image request.
        wcs = 1: only images with world coordinate system (WCS) in their headers
            will be returned. Ignored for a single-image request.
        filter: when present, only images of the specified filter will be
            returned. Ignored for a single-image request.
        explen: when present, only images of the specified exposure length
            (seconds) will be returned. Ignored for a single-image request.
        processed: return specific radio cartographer channel(s)
        channel: if `processed` is set, return the given radio cartographer
            channel ("left", "right", or "composite"); default: return all three
            channels

    :return: for a single-header request, a FITS header in text format; for a
        multi-image request, a pair (filename, data), where `data` is the ZIP
        archive contents
    """
    return api_call(
        'download/header', 'get', kwargs, server=server, api_key=api_key)


def movie(api_key=None, server=None, **kwargs):
    """
    Download images combined into a movie

    :param str api_key: optional API access key
    :param str server: optional API server (for testing)
    :param kwargs::
        image | images: requested image ID(s); multiple IDs are provided as a
            list or a single string, with individual IDs separated by commas.
            An image ID consists of a single-character image type identifier,
            followed by an image ID number. Type Identifiers:
                "r" => raw exposure,
                "m" => master calibration image,
                "w" => Afterglow workspace image,
                "t" => Afterglow temp image,
                "s" => Afterglow sample image
        obs | observations: requested optical observation ID(s); multiple IDs
            are provided as a list or a single string, with individual IDs
            separated by commas
        layer: for multi-HDU observations (e.g. polarimetry or spectral), return
            specific layer (0-based); if omitted, return the primary images
            (layer 0)
        wcs = 1: only images with world coordinate system (WCS) in their headers
            will be returned
        filter: when present, only images of the specified optical filter will
            be returned
        explen: when present, only images of the specified exposure length
            (seconds) will be returned
        reduce: bias, dark, and flat correct the requested image(s)
        reducequiet: same as "reduce" except that in the event of a failure
            during reduction the unreduced image is returned.
        delta: maximum separation (in days) allowed between master calibration
            images and the source image during reduction
        mbias: override the auto-selection of master bias to be used in
            calibration by specifying its ID. Presence of this variable will
            automatically add the option "reduce". If set to 0, this type of
            calibration will not be applied.
        mdark: override the auto-selection of master dark to be used in
            calibration by specifying its ID. Presence of this variable will
            automatically add the option "reduce". If set to 0, this type of
            calibration will not be applied.
        mflat: override the auto-selection of master flat to be used in
            calibration by specifying its ID.  Presence of this variable will
            automatically add the option "reduce". If set to 0, this type of
            calibration will not be applied.
        remove_cosmics = 1: remove cosmic rays from the images; defaults to 1 if
            reduce=1 and to 0 otherwise
        find_cosmics = 1: return images containing only cosmic rays present in
            the original images.
        scale: scale images by the given factor
        width: scale images to the given width
        height: scale images to the given width
        min: specifies the lower percentile of the histogram used to set the
            black point of the image
        max: specifies the upper percentile of the histogram used to set the
            white point of the image
        framerate: specifies the frame rate in frames per second
        keyframeinterval: specifies the key frame interval in seconds
        bitrate: specifies the bitrate in bits per second
        format: specifies the output format

    :return: a pair (filename, data), with `data` containing the movie
    """
    return api_call(
        'download/movie', 'get', kwargs, server=server, api_key=api_key)
