"""
SkyNet API: exceptions
"""


class ApiError(Exception):
    """
    Base class for all SkyNet API exceptions

    Custom exceptions are defined as follows:

        class CustomApiError(ApiError):
            code = 40001
            message = 'Default error message'  # optional

    The unique 5-digit error code consists of a 3-digit HTTP status code (see
    https://en.wikipedia.org/wiki/List_of_HTTP_status_codes) and a 2-digit
    error code for this particular exception.

    An optional default error message may be specified, which allows raising
    exceptions like this:

        raise CustomApiError()

    as well as like this:

        raise CustomApiError('Custom message')
    """
    code = 50000  # full error code = status*100 + subcode
    status = 500  # HTTP status code
    subcode = 0   # unique error subcode for this HTTP status
    message = ''  # default error message

    def __init__(self, message=None):
        self.status, self.code = divmod(self.code, 100)
        if message is not None:
            self.message = message

        super(ApiError, self).__init__(self.message)


class BadRequestError(ApiError):
    code = 40000


class BadAuthenticationDataError(ApiError):
    code = 40001


class InvalidTokenDataError(ApiError):
    code = 40100


class InsufficientPrivilegesError(ApiError):
    code = 40300


class ActionForbiddenError(ApiError):
    code = 40301


class ResourceNotFoundError(ApiError):
    code = 40400


class InternalApiError(ApiError):
    code = 50000


class ReductionFailure(ApiError):
    code = 50001
