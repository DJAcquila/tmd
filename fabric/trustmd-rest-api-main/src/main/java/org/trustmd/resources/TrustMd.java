package org.trustmd.resources;

import org.apache.commons.logging.Log;
import org.apache.commons.logging.LogFactory;
import org.eclipse.microprofile.metrics.annotation.Timed;
import org.eclipse.microprofile.openapi.annotations.OpenAPIDefinition;
import org.eclipse.microprofile.openapi.annotations.Operation;
import org.eclipse.microprofile.openapi.annotations.info.Info;
import org.hyperledger.fabric.gateway.*;
import org.trustmd.model.Trust;

import javax.enterprise.context.ApplicationScoped;
import javax.ws.rs.*;
import javax.ws.rs.core.MediaType;
import javax.ws.rs.core.Response;
import java.io.IOException;
import java.nio.file.Path;
import java.nio.file.Paths;

@javax.ws.rs.Path("/trustmd")
@ApplicationScoped
@OpenAPIDefinition(info = @Info(title = "trustmd endpoint", version = "1.0"))
public class TrustMd {
    private static final Log _logger = LogFactory.getLog(TrustMd.class);
    private final static Path CURRENT_HOME_PATH = Paths.get(System.getProperty("user.home"));
    private final static Path ROOT_PATH = Paths.get(CURRENT_HOME_PATH.toString(), "workspace", "minifab", "trustmdtest");
    private final static Path CONNECTIONS_PATH = Paths.get(ROOT_PATH.toString(), "Connections");
    private final static String CONNECTIONS_FILE = "/v1-3p-connection.json";

    static {
        System.setProperty("org.hyperledger.fabric.sdk.service_discovery.as_localhost", "true");
    }

    @Timed(name = "createTrustTime", tags = {
            "method=post" }, absolute = true, description = "Time needed to add a trust value to the ledger")
    @POST
    @Produces(MediaType.APPLICATION_JSON)
    @Consumes(MediaType.APPLICATION_JSON)
    @javax.ws.rs.Path("/createTrust")
    @Operation(summary = "Add a new Trust info to the ledger", description = "Requires a unique key to be successfully executed")
    public Response createTrust(Trust aTrust) {
        try {
            System.out.println("Calling CreateTrust");
            final Gateway.Builder gBuilder = getGatewayBuilder();

            try (final Gateway gateway = gBuilder.connect()) {
                final Contract contract = getContract(gateway);
                final byte[] result = contract.submitTransaction("createTrust",
                        aTrust.getClusterHeadId(),
                        aTrust.getEvaluatedNodeId(),
                        aTrust.getTrustDecision(),
                        aTrust.getDegreeOfBelief(),
                        aTrust.getMecHostId());
                return Response.status(Response.Status.OK).entity(new String(result)).build();
            } catch (Exception exception) {
                throw new javax.ws.rs.ServerErrorException("Unable to get network/contract and execute query", Response.Status.INTERNAL_SERVER_ERROR);
            }
        } catch (Exception exception) {
            // Initially raise ServiceUnavailable
            throw new javax.ws.rs.ServerErrorException("Unable to find config or wallet", Response.Status.NOT_FOUND);
        }
    }

    @Timed(name = "updateTrustTime", tags = {
            "method=put" }, absolute = true, description = "Time needed to update a trust value to the ledger")
    @PUT
    @Produces(MediaType.APPLICATION_JSON)
    @Consumes(MediaType.APPLICATION_JSON)
    @javax.ws.rs.Path("/update")
    @Operation(summary = "Update a new Trust info to the ledger", description = "Requires a unique key to be successfully executed")
    public Response updateTrust(final Trust aTrust) {
        try {
            System.out.println("Calling UpdateTrust");
            final Gateway.Builder gBuilder = getGatewayBuilder();

            try (final Gateway gateway = gBuilder.connect()) {
                final Contract contract = getContract(gateway);
                final byte[] result = contract.submitTransaction("updateTrust",
                        aTrust.getClusterHeadId(),
                        aTrust.getEvaluatedNodeId(),
                        aTrust.getTrustDecision(),
                        aTrust.getDegreeOfBelief(),
                        aTrust.getMecHostId());
                return Response.status(Response.Status.OK).entity(new String(result)).build();
            } catch (Exception exception) {
                throw new javax.ws.rs.ServerErrorException("Unable to get network/contract and execute query", Response.Status.INTERNAL_SERVER_ERROR);
            }
        } catch (Exception exception) {
            // Initially raise ServiceUnavailable
            throw new javax.ws.rs.ServerErrorException("Unable to find config or wallet", Response.Status.NOT_FOUND);
        }
    }

    @Timed(name = "getTrustTime", tags = {
            "method=get" }, absolute = true, description = "Time needed to get a trust value of the ledger")
    @GET
    @Produces(MediaType.APPLICATION_JSON)
    @Consumes(MediaType.APPLICATION_JSON)
    @javax.ws.rs.Path("{id}")
    @Operation(summary = "Get the Trust info from the ledger", description = "Requires a unique key to be successfully executed")
    public Response getTrustByNodeId(@PathParam("id") String evaluatedNodeId) {
        try {
            System.out.println("Calling GetTrust " + evaluatedNodeId);
            final Gateway.Builder gBuilder = getGatewayBuilder();

            try (final Gateway gateway = gBuilder.connect()) {
                final Contract contract = getContract(gateway);
                final byte[] result = contract.submitTransaction("getTrust", evaluatedNodeId);

                return Response.status(Response.Status.OK).entity(new String(result)).build();
            } catch (Exception exception) {
                throw new javax.ws.rs.ServerErrorException("Unable to get network/contract and execute query", Response.Status.INTERNAL_SERVER_ERROR);
            }
        } catch (Exception exception) {
            // Initially raise ServiceUnavailable
            throw new javax.ws.rs.ServerErrorException("Unable to find config or wallet", Response.Status.NOT_FOUND);
        }
    }

    @Timed(name = "deleteTrustTime", tags = {
            "method=delete" }, absolute = true, description = "Time needed to get a trust value of the ledger")
    @DELETE
    @Produces(MediaType.APPLICATION_JSON)
    @Consumes(MediaType.APPLICATION_JSON)
    @javax.ws.rs.Path("{id}")
    @Operation(summary = "Delete Trust info from the ledger", description = "Requires a unique key to be successfully executed")
    public Response deleteByNodeId(@PathParam("id") String evaluatedNodeId) {
        try {
            System.out.println("Calling DeleteTrust " + evaluatedNodeId);
            final Gateway.Builder gBuilder = getGatewayBuilder();

            try (final Gateway gateway = gBuilder.connect()) {
                final Contract contract = getContract(gateway);
                final byte[] result = contract.submitTransaction("delete", evaluatedNodeId);

                return Response.status(Response.Status.OK).entity(new String(result)).build();
            } catch (Exception exception) {
                throw new javax.ws.rs.ServerErrorException("Unable to get network/contract and execute query", Response.Status.INTERNAL_SERVER_ERROR);
            }
        } catch (Exception exception) {
            // Initially raise ServiceUnavailable
            throw new javax.ws.rs.ServerErrorException("Unable to find config or wallet", Response.Status.NOT_FOUND);
        }
    }
    @GET
    @Produces(MediaType.APPLICATION_JSON)
    @Consumes(MediaType.APPLICATION_JSON)
    @javax.ws.rs.Path("/walletTest")
    public Response walletTest(Trust aTrust) throws IOException {
        // expecting wallet directory within the default server location
        // wallet exported from Fabric wallets
        final Gateway.Builder gBuilder = getGatewayBuilder();
        try (final Gateway gateway = gBuilder.connect()) {
            final Contract contract = getContract(gateway);
            System.out.println(contract);
        }

        return Response.status(Response.Status.OK).entity("OK").build();
    }

    private Wallet getWallet() throws IOException {
        final Path walletPath = Paths.get(ROOT_PATH.toString(), "vars", "app", "go", "wallets");
        return Wallets.newFileSystemWallet(walletPath);
    }

    private Gateway.Builder getGatewayBuilder() throws IOException {
        final Wallet wallet = getWallet();
        final Identity adminIdentity = wallet.get("Admin");
        if (adminIdentity != null) {
            System.out.println("Entity adminIdentity already registered");
        }

        // Loading connections config file from expected path
        final Path networkConfigPath = Paths.get(CONNECTIONS_PATH + CONNECTIONS_FILE);
        System.out.println("Network config path: " + networkConfigPath);

        return Gateway.createBuilder()
                .identity(wallet, "Admin")
                .networkConfig(networkConfigPath)
                .discovery(false);
    }

    private Contract getContract(final Gateway gateway) {
        final Network network = gateway.getNetwork("mychannel");
        return network.getContract("trustmd");
    }
}
