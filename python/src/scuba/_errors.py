class ScubaError(Exception):
    """Raised when SCUBA.jl raises one of its typed errors (VarsNotFound, MetadataNotFound, ...)."""

    def __init__(self, message: str, kind: str | None = None):
        super().__init__(message)
        self.kind = kind


class ScubaWarning(UserWarning):
    """A structured SCUBA warning re-emitted from Julia; `.code` holds the Julia warning code."""

    def __init__(self, message: str, code: str | None = None, vars: list[str] | None = None):
        super().__init__(message)
        self.code = code
        self.vars = vars or []
