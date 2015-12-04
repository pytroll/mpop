Feature: Simple and intuitive scene loading

    The scientific user explores the data and prototypes new algorithms. 
    It needs access not only to the calibrated data,
    but also to the raw data and probably a majority of the metadata. 
    The user would work with data locally, and it has to be easy to tell mpop2 where the data is. 
    Providing filename templates or editing config file before starting working is a pain, so it should be avoided.

    To load the data should be a simple 1-step procedure. 
    At load time, the user provides the data and metadata he/she needs, 
    and if some items are unavailable/unaccessible, 
    the user should be informed in a gentle but clear way (ie. no crash). 
    The data and metadata available from the file have to be explorable, 
    so that the user don’t need to guess what the (meta)data is called.

    Scenario: 1-step data loading
        Given data is available
        When user loads the data without providing a config file
        Then scene is returned

    Scenario: No crash when metadata is missing
        Given data is available
        When user loads the data without providing a config file
        And some items are no available
        Then scene is returned
