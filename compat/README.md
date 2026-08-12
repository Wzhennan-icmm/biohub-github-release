# Legacy compatibility wrappers

These v1 wrappers preserve common historic Python entry names while dispatching
to BioHub. They are deprecated: use `biohub <group> <command>` or
`biohub run <script-id>` in new workflows. Keep wrappers through v1.x and cover
them with regression fixtures before any v2 removal.
